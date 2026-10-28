# Nonlinear operator updates in the original full/global field space.
_has_transient_modifiers(callbacks...) = any(f -> f !== nothing, callbacks)

function _transient_nonlinear_options(options)
    options isa NamedTuple || throw(ArgumentError("nonlinearoptions must be a NamedTuple."))
    all(k -> k in (:maxiter,:reltol,:abstol,:relaxation), keys(options)) ||
        throw(ArgumentError("Unknown nonlinear option; use maxiter, reltol, abstol, relaxation."))
    maxiter=get(options,:maxiter,30)
    reltol=get(options,:reltol,1e-9)
    abstol=get(options,:abstol,1e-11)
    relaxation=get(options,:relaxation,1.0)
    maxiter isa Integer && !(maxiter isa Bool) && maxiter>0 ||
        throw(ArgumentError("nonlinear maxiter must be a positive integer."))
    all(x -> x isa Real && isfinite(x), (reltol,abstol,relaxation)) &&
        reltol>=0 && abstol>=0 && 0<relaxation<=1 ||
        throw(ArgumentError("Invalid nonlinear tolerances/relaxation."))
    return (;maxiter,reltol,abstol,relaxation)
end

function _transient_operator_copy(A)
    A === nothing && return nothing
    P=parent_system(A)
    B=SystemMatrix(copy(P.A),P.model,P.test_model,P.problems,P.offsets)
    return A isa SymmetricSystemMatrix ? Symmetric(B,A.uplo) : B
end

function _transient_zero_operator(p)
    K=p.layout
    return SystemMatrix(spzeros(size(K.A)...),K.model,K.test_model,K.problems,K.offsets)
end

function _transient_rhs_field(values,p,t)
    fields=split_multifield_solution(reshape(values,:,1),p.problems,p.offsets,[t])
    return p.layout.problems === nothing ? fields :
        SystemVector(reshape(Vector(values),:,1),nothing,p.problems,p.offsets)
end

function _transient_snapshot(p,x,v,a,j)
    function fields(free_values,boundary)
        z=zeros(size(p.T,2),1)
        z[p.free,1].=free_values
        z[p.fixed,1].=boundary
        result=split_multifield_solution(p.L*z,p.problems,p.offsets,[p.t[j]])
        return result isa Tuple ? result : (result,)
    end
    return (;fields=fields(x,p.D[:,j]),
            velocity=v === nothing ? nothing : fields(v[1],v[2]),
            acceleration=a === nothing ? nothing : fields(a[1],a[2]),
            problems=Tuple(p.problems),t=p.t[j],step=j-1)
end

function _transient_initial_frame(K,M,C,p)
    op=(A2=M,A1=C === nothing ? Symmetric(_transient_zero_operator(p)) : C,A0=K)
    return (;operators=op,A2=p.M,A1=p.C,A0=p.K,
            b=p.T'* (p.Q === nothing ? p.Ffull[:,1] : p.Q'*p.Ffull[:,1]))
end

function _transient_same_operator(A,B)
    A === nothing && return B === nothing
    B === nothing && return false
    return (A isa SymmetricSystemMatrix)==(B isa SymmetricSystemMatrix) &&
        (!(A isa SymmetricSystemMatrix) || A.uplo==B.uplo) &&
        parent_system(A).A == parent_system(B).A
end

function _transient_modified_frame(base,latest,p,snapshot,previous,j,iteration,phase,hooks,dt)
    operators=base.operators
    matrices=(A2=latest.A2,A1=latest.A1,A0=latest.A0)
    # New inertia/stiffness are visible to damping (e.g. Rayleigh damping).
    for key in (:A2,:A0,:A1)
        callback=getproperty(hooks,key)
        callback === nothing && continue
        state=merge(snapshot,(;previous,iteration,phase,dt,operators))
        argument=_transient_operator_copy(getproperty(base.operators,key))
        applicable(callback,argument,state) ||
            throw(ArgumentError("modify$key must be callable as (operator, state)."))
        output=callback(argument,state)
        value=output === nothing ? getproperty(base.operators,key) : output
        value isa _TransientMatrix || throw(ArgumentError("modify$key must return a SystemMatrix, Symmetric(SystemMatrix), or nothing."))
        _transient_check_matrix(p.layout,parent_system(value))
        old=getproperty(latest.operators,key)
        projected=_transient_same_operator(value,old) ? getproperty(latest,key) :
            p.L'*(parent_system(value).A*p.L)
        operators=merge(operators,NamedTuple{(key,)}((value,)))
        matrices=merge(matrices,NamedTuple{(key,)}((projected,)))
    end
    global_rhs=Vector(_transient_column(p.Ffull,j))
    if hooks.RHS !== nothing
        state=merge(snapshot,(;previous,iteration,phase,dt,operators))
        rhs=_transient_rhs_field(global_rhs,p,p.t[j])
        applicable(hooks.RHS,rhs,state) ||
            throw(ArgumentError("modifyRHS must be callable as (rhs, state)."))
        output=hooks.RHS(rhs,state)
        if output !== nothing
            # All callback fields, matrices and RHS components are global.
            rhs_values=_transient_values(output,p.layout,nothing,1)
            global_rhs=Vector(@view(rhs_values[:,1]))
        end
    end
    b=p.L'*global_rhs
    inactive=p.inactive
    if !isempty(inactive)
        any(!iszero,b[inactive]) && throw(ArgumentError("modifyRHS activates an inactive MPC DOF."))
        for A in values(matrices)
            (nnz(A[inactive,:])==0 || all(iszero,nonzeros(A[inactive,:]))) &&
            (nnz(A[:,inactive])==0 || all(iszero,nonzeros(A[:,inactive]))) ||
                throw(ArgumentError("modify callback activates an inactive MPC DOF; the active layout must stay fixed."))
        end
    end
    return merge((;operators,b),matrices)
end

# Matrix values and symmetry declarations form the cache key. Returning a new
# but equal matrix (or changing only RHS) does not rebuild a factor/preconditioner.
function _transient_modified_prepare!(cache,key,A,operators,solveroptions;mass=false,cdm=false)
    symmetric=all(A -> A === nothing || A isa SymmetricSystemMatrix,operators)
    diagonal=mass && isdiag(A)
    entry=get(cache,key,nothing)
    if entry === nothing || entry.A!=A || entry.symmetric!=symmetric
        factor=(cdm || diagonal) ? _transient_factor(A;mass) : nothing
        prepared=factor === nothing ? _prepare_transient_linear_solver(A,operators;solveroptions...) : nothing
        entry=(;A=copy(A),symmetric,factor,prepared)
        cache[key]=entry
    end
    return entry
end

_transient_modified_solve(entry,b) = entry.factor === nothing ?
    solve_prepared_linear_system(entry.prepared,Vector(b)) : entry.factor\Vector(b)

function _transient_modified_solve!(cache,key,A,operators,b,solveroptions;kwargs...)
    entry=_transient_modified_prepare!(cache,key,A,operators,solveroptions;kwargs...)
    return _transient_modified_solve(entry,b)
end

function _transient_nonlinear_converged(terms,options;reference=0.0)
    residual=sum(terms)
    all(isfinite,residual) || error("Transient nonlinear residual is not finite.")
    scale=sum(norm,terms)+reference
    return norm(residual)<=options.abstol+options.reltol*scale,norm(residual)
end

function _transient_nonlinear_failure(method,p,j,residual,options)
    error("$method: nonlinear iteration did not converge at step $(j-1), t=$(p.t[j]); " *
          "residual=$residual after $(options.maxiter) updates.")
end

# At fixed displacement/velocity, determine acceleration, including any user
# acceleration dependence. CDM remains explicit with respect to displacement.
function _transient_modified_acceleration(p,base,x,v,guess,j,Vd,Ad,previous,hooks,options,cache,linear;
                                           cdm=false,phase=:initial,dt)
    a=copy(guess); frame=base
    for iteration in 0:options.maxiter
        snapshot=_transient_snapshot(p,x,(v,Vd[:,j]),(a,Ad[:,j]),j)
        frame=_transient_modified_frame(base,frame,p,snapshot,previous,j,iteration,phase,hooks,dt)
        f,c=p.free,p.fixed
        rhs=frame.b[f]-frame.A0[f,f]*x-frame.A0[f,c]*p.D[:,j]-
            frame.A1[f,f]*v-frame.A1[f,c]*Vd[:,j]-frame.A2[f,c]*Ad[:,j]
        # Validate/prepare the free mass even for a zero initial residual.
        candidate=isempty(f) ? Float64[] : _transient_modified_solve!(cache,:mass,
            frame.A2[f,f],(frame.operators.A2,),rhs,linear;mass=true,cdm)
        reference=norm(frame.b[f])+norm(frame.A0[f,f]*x)+norm(frame.A0[f,c]*p.D[:,j])+
            norm(frame.A1[f,f]*v)+norm(frame.A1[f,c]*Vd[:,j])+norm(frame.A2[f,c]*Ad[:,j])
        ok,residual=_transient_nonlinear_converged((frame.A2[f,f]*a,-rhs),options;reference)
        ok && return a,frame,snapshot
        iteration==options.maxiter && _transient_nonlinear_failure(cdm ? "CDM" : "HHT",p,j,residual,options)
        a .+= options.relaxation.*(candidate-a)
    end
end

function _modified_hht(K,M,C,rhs,U0,V0,n,dt;α,γ,β,solver,preconditioner,
    solveroptions,preconditioneroptions,modifyA2,modifyA1,modifyA0,modifyRHS,
    nonlinearoptions,kwargs...)
    options=_transient_nonlinear_options(nonlinearoptions)
    p=prepare_transient_system(K,M,C,rhs,U0,V0,n,dt;_modifier_data=true,kwargs...)
    hooks=(A2=modifyA2,A1=modifyA1,A0=modifyA0,RHS=modifyRHS)
    linear=(;solver,preconditioner,solveroptions,preconditioneroptions)
    cache=Dict{Symbol,Any}()
    Vd,Ad=_transient_derivatives(p.D,dt)
    for j in 2:n
        Ad[:,j].=(p.D[:,j]-p.D[:,j-1]-dt*Vd[:,j-1]-dt^2*(0.5-β)*Ad[:,j-1])/(β*dt^2)
        Vd[:,j].=Vd[:,j-1]+dt*((1-γ)*Ad[:,j-1]+γ*Ad[:,j])
    end
    U,V=zeros(length(p.free),n),zeros(length(p.free),n)
    U[:,1],V[:,1]=p.u0,p.v0
    base=_transient_initial_frame(K,M,C,p)
    a,base,previous=_transient_modified_acceleration(p,base,p.u0,p.v0,zeros(length(p.free)),1,
        Vd,Ad,nothing,hooks,options,cache,linear;dt)
    f,c=p.free,p.fixed
    for j in 2:n
        xold,vold=copy(U[:,j-1]),copy(V[:,j-1])
        oldforce=base.A0[f,f]*xold+base.A0[f,c]*p.D[:,j-1]+
                 base.A1[f,f]*vold+base.A1[f,c]*Vd[:,j-1]-base.b[f]
        oldreference=norm(base.b[f])+norm(base.A0[f,c]*p.D[:,j-1])+norm(base.A1[f,c]*Vd[:,j-1])
        xp=xold+dt*vold+dt^2*(0.5-β)*a
        vp=vold+dt*(1-γ)*a
        trial=copy(a); frame=base
        for iteration in 0:options.maxiter
            x=xp+β*dt^2*trial; v=vp+γ*dt*trial
            snapshot=_transient_snapshot(p,x,(v,Vd[:,j]),(trial,Ad[:,j]),j)
            frame=_transient_modified_frame(base,frame,p,snapshot,previous,j,iteration,:step,hooks,dt)
            inertia=frame.A2[f,f]*trial+frame.A2[f,c]*Ad[:,j]
            force=frame.A0[f,f]*x+frame.A0[f,c]*p.D[:,j]+
                  frame.A1[f,f]*v+frame.A1[f,c]*Vd[:,j]-frame.b[f]
            terms=(inertia,(1+α)*force,-α*oldforce)
            A=frame.A2[f,f]+(1+α)*γ*dt*frame.A1[f,f]+(1+α)*β*dt^2*frame.A0[f,f]
            entry=isempty(f) ? nothing : _transient_modified_prepare!(cache,:effective,A,
                values(frame.operators),linear)
            reference=(1+α)*(norm(frame.b[f])+norm(frame.A0[f,c]*p.D[:,j])+
                norm(frame.A1[f,c]*Vd[:,j]))+norm(frame.A2[f,c]*Ad[:,j])+abs(α)*oldreference
            ok,residual=_transient_nonlinear_converged(terms,options;reference)
            if ok
                U[:,j],V[:,j]=x,v
                a=trial; base=frame; previous=snapshot
                break
            end
            iteration==options.maxiter && _transient_nonlinear_failure("HHT",p,j,residual,options)
            correction=_transient_modified_solve(entry,-sum(terms))
            trial .+= options.relaxation.*correction
        end
    end
    return reconstruct_transient_solution(U,p.D,p),reconstruct_transient_solution(V,Vd,p)
end

function _modified_cdm(K,M,C,rhs,U0,V0,n,dt;modifyA2,modifyA1,modifyA0,modifyRHS,
                       nonlinearoptions,kwargs...)
    options=_transient_nonlinear_options(nonlinearoptions)
    p=prepare_transient_system(K,M,C,rhs,U0,V0,n,dt;_modifier_data=true,kwargs...)
    hooks=(A2=modifyA2,A1=modifyA1,A0=modifyA0,RHS=modifyRHS)
    cache=Dict{Symbol,Any}(); Vd,Ad=_transient_derivatives(p.D,dt)
    U,V=zeros(length(p.free),n),zeros(length(p.free),n)
    U[:,1],V[:,1]=p.u0,p.v0
    base=_transient_initial_frame(K,M,C,p)
    a,base,previous=_transient_modified_acceleration(p,base,p.u0,p.v0,zeros(length(p.free)),1,
        Vd,Ad,nothing,hooks,options,cache,(;);cdm=true,dt)
    prior=p.u0-dt*p.v0+dt^2/2*a
    for j in 1:n-1
        x=copy(U[:,j]); v=(x-prior)/dt
        a,base,previous=_transient_modified_acceleration(p,base,x,v,a,j,Vd,Ad,previous,
            hooks,options,cache,(;);cdm=true,phase=:explicit,dt)
        U[:,j+1].=2x-prior+dt^2*a
        V[:,j+1].=(U[:,j+1]-x)/dt
        prior=x
    end
    return reconstruct_transient_solution(U,p.D,p),reconstruct_transient_solution(V,Vd,p)
end

function _modified_fdm(K,C,rhs,U0,n,dt;ϑ,solver,preconditioner,solveroptions,
    preconditioneroptions,modifyA1,modifyA0,modifyRHS,nonlinearoptions,kwargs...)
    options=_transient_nonlinear_options(nonlinearoptions)
    p=prepare_transient_system(K,nothing,C,rhs,U0,nothing,n,dt;_modifier_data=true,kwargs...)
    hooks=(A2=nothing,A1=modifyA1,A0=modifyA0,RHS=modifyRHS)
    linear=(;solver,preconditioner,solveroptions,preconditioneroptions)
    cache=Dict{Symbol,Any}(); U=zeros(length(p.free),n); U[:,1]=p.u0
    base=_transient_initial_frame(K,nothing,C,p)
    previous=_transient_snapshot(p,p.u0,nothing,nothing,1)
    base=_transient_modified_frame(base,base,p,previous,nothing,1,0,:initial,hooks,dt)
    f,c=p.free,p.fixed
    for j in 2:n
        old=copy(U[:,j-1]); trial=copy(old); frame=base
        oldforce=base.A0[f,f]*old+base.A0[f,c]*p.D[:,j-1]-base.b[f]
        for iteration in 0:options.maxiter
            rate=(trial-old)/dt; boundary_rate=(p.D[:,j]-p.D[:,j-1])/dt
            snapshot=_transient_snapshot(p,trial,(rate,boundary_rate),nothing,j)
            frame=_transient_modified_frame(base,frame,p,snapshot,previous,j,iteration,:step,hooks,dt)
            capacity=(1-ϑ)*base.A1+ϑ*frame.A1
            force=frame.A0[f,f]*trial+frame.A0[f,c]*p.D[:,j]-frame.b[f]
            terms=(capacity[f,f]*rate+capacity[f,c]*boundary_rate,ϑ*force,(1-ϑ)*oldforce)
            A=capacity[f,f]+ϑ*dt*frame.A0[f,f]
            entry=isempty(f) ? nothing : _transient_modified_prepare!(cache,:effective,A,
                (base.operators.A1,frame.operators.A1,frame.operators.A0),linear)
            reference=ϑ*(norm(frame.b[f])+norm(frame.A0[f,c]*p.D[:,j]))+
                (1-ϑ)*(norm(base.b[f])+norm(base.A0[f,c]*p.D[:,j-1]))+
                norm(capacity[f,c]*boundary_rate)
            ok,residual=_transient_nonlinear_converged(terms,options;reference)
            if ok
                U[:,j]=trial; base=frame; previous=snapshot
                break
            end
            iteration==options.maxiter && _transient_nonlinear_failure("FDM",p,j,residual,options)
            correction=_transient_modified_solve(entry,-dt*sum(terms))
            trial .+= options.relaxation.*correction
        end
    end
    return reconstruct_transient_solution(U,p.D,p)
end
