@testset "modify operators: linear equivalence and callback state" begin
    with_test_gmsh() do
        line_mesh(n=2,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        M=∫(P⋅P); K=∫(Grad(P)⋅Grad(P))+M; C=0.2M
        n,dt=17,0.01; t=collect(0:n-1)*dt; N=ndofs(P)
        zero=transient_field(P,zeros(N,1),[0.0])
        load=transient_field(P,repeat(reshape(sin.(t),1,:),N,1),t)
        support=[BoundaryCondition("left";field=P,p=(x,y,z,t)->0.01t)]
        seen=Any[]
        function unchanged(A,s)
            push!(seen,(s.step,s.iteration,s.phase,s.t))
            @test s.fields isa Tuple && length(s.fields)==1
            @test s.fields[1].model === s.problems[1] === P
            @test s.fields[1].nsteps==1 && s.fields[1].t==[s.t]
            @test s.dt ≈ dt
            return nothing
        end
        for run in (HHT,CDM)
            u,v=run(K,M,C,load,zero,zero,n,dt;support)
            um,vm=run(K,M,C,load,zero,zero,n,dt;support,
                modifyA2=unchanged,modifyA1=(A,s)->nothing,
                modifyA0=(A,s)->nothing,modifyRHS=(b,s)->nothing)
            @test um.a ≈ u.a atol=1e-10
            @test vm.a ≈ v.a atol=1e-9
        end
        x=FDM(K,M,load,zero,n,dt;support)
        xm=FDM(K,M,load,zero,n,dt;support,modifyA1=unchanged)
        @test xm.a ≈ x.a atol=1e-10
        @test any(s->s[2]>0,seen)
        @test any(s->s[1]==0,seen)
        @test_throws ArgumentError FDM(K,M,load,zero,n,dt;modifyA2=unchanged)
    end
end

@testset "modify operators: time-dependent manufactured solutions" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        N=ndofs(P); I0=SystemMatrix(spdiagm(0=>ones(N)),P)
        matrix(value)=SystemMatrix(spdiagm(0=>fill(value,N)),P)
        field(value,t)=transient_field(P,fill(value,N,1),[t])
        n,dt=31,0.01; t=collect(0:n-1)*dt
        a2(t)=1+0.1t; a1(t)=0.2+0.05t; a0(t)=2+0.3t
        u(t)=0.5+0.3t+t^2; v(t)=0.3+2t
        mass=(A,s)->matrix(a2(s.t))
        stiffness=(A,s)->matrix(a0(s.t))
        damping=(A,s)->begin
            @test LowLevelFEM.parent_system(s.operators.A2).A == matrix(a2(s.t)).A
            @test LowLevelFEM.parent_system(s.operators.A0).A == matrix(a0(s.t)).A
            matrix(a1(s.t))
        end
        force=(b,s)->field(a2(s.t)*2+a1(s.t)*v(s.t)+a0(s.t)*u(s.t),s.t)
        zero=field(0.0,0.0); initial=field(u(0.0),0.0); velocity=field(v(0.0),0.0)
        x,xd=HHT(2I0,I0,0.2I0,zero,initial,velocity,n,dt;
            modifyA2=mass,modifyA1=damping,modifyA0=stiffness,modifyRHS=force,copyOperators=false)
        @test x.a ≈ repeat(reshape(u.(t),1,:),N,1) atol=2e-10
        @test xd.a ≈ repeat(reshape(v.(t),1,:),N,1) atol=2e-9
        explicit_force=(b,s)->field(a2(s.t)*2+a1(s.t)*(v(s.t)-(s.phase==:initial ? 0 : dt))+
                                    a0(s.t)*u(s.t),s.t)
        x,_=CDM(2I0,I0,0.2I0,zero,initial,velocity,n,dt;
            modifyA2=mass,modifyA1=damping,modifyA0=stiffness,modifyRHS=explicit_force,copyOperators=false)
        @test x.a ≈ repeat(reshape(u.(t),1,:),N,1) atol=2e-10
        for θ in (0.0,0.5,1.0)
            x=FDM(2I0,0.2I0,zero,initial,n,dt;ϑ=θ,
                modifyA1=(A,s)->matrix(a1(s.t)),modifyA0=stiffness,
                modifyRHS=(b,s)->field(a1(s.t)+a0(s.t)*(0.5+s.t),s.t),copyOperators=false)
            @test x.a ≈ repeat(reshape(0.5 .+t,1,:),N,1) atol=2e-10
        end
    end
end

@testset "modify operators: nonlinear HHT/FDM/CDM" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        N=ndofs(P); I0=SystemMatrix(spdiagm(0=>ones(N)),P)
        zero=transient_field(P,zeros(N,1),[0.0]); initial=transient_field(P,ones(N,1),[0.0])
        k,λ,c=2.0,4.0,0.1
        n,dt=41,0.025; iterations=Int[]
        function tangent(A,s)
            push!(iterations,s.iteration)
            values=vec(s.fields[1].a)
            return SystemMatrix(spdiagm(0=>k .+3λ.*values.^2),P)
        end
        forcing=(b,s)->transient_field(P,2λ.*s.fields[1].a.^3,[s.t])
        U,V=HHT(k*I0,I0,c*I0,zero,initial,zero,n,dt;
            modifyA0=tangent,modifyRHS=forcing,
            nonlinearoptions=(reltol=1e-12,abstol=1e-13,))
        borrowedU,borrowedV=HHT(k*I0,I0,c*I0,zero,initial,zero,n,dt;
            modifyA0=tangent,modifyRHS=forcing,copyOperators=false,
            nonlinearoptions=(reltol=1e-12,abstol=1e-13,))
        @test borrowedU.a ≈ U.a atol=1e-12
        @test borrowedV.a ≈ V.a atol=1e-12
        x,v=1.0,0.0; a=-k-λ
        for j in 2:n
            xp=x+dt*v+dt^2/4*a; vp=v+dt/2*a
            trial=a
            for it in 1:15
                xx=xp+dt^2/4*trial; vv=vp+dt/2*trial
                r=trial+c*vv+k*xx+λ*xx^3
                abs(r)<1e-13 && break
                trial-=r/(1+c*dt/2+(k+3λ*xx^2)*dt^2/4)
            end
            a=trial; x=xp+dt^2/4*a; v=vp+dt/2*a
            @test U.a[:,j] ≈ fill(x,N) atol=1e-11
            @test V.a[:,j] ≈ fill(v,N) atol=1e-10
        end
        @test maximum(iterations)>=2
        X=FDM(k*I0,I0,zero,initial,n,dt;ϑ=1.0,
            modifyA0=tangent,modifyRHS=forcing,
            nonlinearoptions=(reltol=1e-12,abstol=1e-13,))
        borrowedX=FDM(k*I0,I0,zero,initial,n,dt;ϑ=1.0,
            modifyA0=tangent,modifyRHS=forcing,copyOperators=false,
            nonlinearoptions=(reltol=1e-12,abstol=1e-13,))
        @test borrowedX.a ≈ X.a atol=1e-12
        x=1.0
        for j in 2:n
            trial=x
            for it in 1:15
                r=(trial-x)/dt+k*trial+λ*trial^3
                abs(r)<1e-13 && break
                trial-=r/(1/dt+k+3λ*trial^2)
            end
            x=trial
            @test X.a[:,j] ≈ fill(x,N) atol=1e-11
        end
        U,V=CDM(k*I0,I0,c*I0,zero,initial,zero,n,dt;
            modifyA0=tangent,modifyRHS=forcing)
        borrowedU,borrowedV=CDM(k*I0,I0,c*I0,zero,initial,zero,n,dt;
            modifyA0=tangent,modifyRHS=forcing,copyOperators=false)
        @test borrowedU.a ≈ U.a atol=1e-12
        @test borrowedV.a ≈ V.a atol=1e-12
        x=1.0; prior=x+dt^2/2*(-k-λ)
        for j in 2:n
            xn=2x-prior+dt^2*(-k*x-λ*x^3-c*(x-prior)/dt)
            @test U.a[:,j] ≈ fill(xn,N) atol=1e-10
            prior,x=x,xn
        end
        saved=copy(I0.A),copy(initial.a)
        @test_throws ErrorException HHT(k*I0,I0,c*I0,zero,initial,zero,3,0.5;
            modifyA0=tangent,modifyRHS=forcing,nonlinearoptions=(maxiter=1,abstol=1e-14,reltol=1e-14))
        @test saved==(I0.A,initial.a)
        @test_throws ArgumentError HHT(I0,I0,zero,initial,zero,n,dt;modifyA0=4)
        @test_throws ArgumentError HHT(I0,I0,zero,initial,zero,n,dt;
            modifyA0=(A,s)->A.A)
        @test_throws ArgumentError HHT(I0,I0,zero,initial,zero,n,dt;
            modifyA0=tangent,nonlinearoptions=(maxiter=0,))
    end
end

@testset "modify operators: remote multifield and inactive validation" begin
    with_test_gmsh() do
        f=transient_remote_fixture(2)
        (;U,Φ,K,M,C,initial,support,mpc,dt,R)=f
        n=7
        force=SystemVector([∫(U⋅[0.0,1.0];Γ="right"),∫(Φ⋅0.0;Ω="body")])
        exact,_=HHT(2K,M,C,force,initial,initial,n,dt;support,mpc)
        seen=Ref(false)
        callback=(A,s)->begin
            @test s.problems==(U,Φ)
            @test length(s.fields)==2
            seen[]=true
            2K
        end
        trial,_=HHT(K,M,C,force,initial,initial,n,dt;support,mpc,modifyA0=callback)
        borrowed,_=HHT(K,M,C,force,initial,initial,n,dt;
            support,mpc,modifyA0=callback,copyOperators=false)
        @test trial[1].a+R.A*trial[2].a ≈ exact[1].a+R.A*exact[2].a atol=1e-10
        @test borrowed[1].a+R.A*borrowed[2].a ≈ exact[1].a+R.A*exact[2].a atol=1e-10
        @test seen[]
        bad=(A,s)->SystemMatrix(LowLevelFEM.parent_system(M).A+spdiagm(0=>ones(size(M.A,1))),
                               M.model,M.test_model,M.problems,M.offsets)
        @test_throws ArgumentError HHT(K,M,C,force,initial,initial,n,dt;
            support,mpc,modifyA2=bad)
    end
end

@testset "modify operators: unilateral contact law" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        N=ndofs(P); I0=SystemMatrix(spdiagm(0=>ones(N)),P)
        zero=transient_field(P,zeros(N,1),[0.0]); speed=transient_field(P,ones(N,1),[0.0])
        k,kc,gap,c=1.0,200.0,0.1,0.05; n,dt=121,0.005
        active=Bool[]
        kt=(A,s)->begin
            mask=vec(s.fields[1].a).>gap
            append!(active,mask)
            SystemMatrix(spdiagm(0=>k .+kc.*mask),P)
        end
        load=(b,s)->transient_field(P,kc*gap.*(s.fields[1].a.>gap),[s.t])
        u,v=HHT(k*I0,I0,c*I0,zero,zero,speed,n,dt;
            modifyA0=kt,modifyRHS=load,nonlinearoptions=(reltol=1e-12,abstol=1e-13))
        x,vel,a=0.0,1.0,-c
        for j in 2:n
            xp=x+dt*vel+dt^2/4*a; vp=vel+dt/2*a
            base=4/dt^2+2c/dt
            numerator=base*xp-c*vp
            trial=numerator/(base+k)
            x=trial<=gap ? trial : (numerator+kc*gap)/(base+k+kc)
            a=(x-xp)*4/dt^2; vel=vp+dt/2*a
            @test u.a[:,j] ≈ fill(x,N) atol=1e-10
            @test v.a[:,j] ≈ fill(vel,N) atol=1e-9
        end
        @test any(active) && any(!,active)
        @test minimum(v.a)<0
        @test u.a[1,end]<gap
    end
end

@testset "modify operators: global coordinates and reduced/spectral projection" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1,order=1)
        P=Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        N=ndofs(P); M=SystemMatrix(spdiagm(0=>ones(N)),P); K=2M; C=0.1M
        cs=[CoordinateSystem(VectorField(P,"body",[0.0,1.0,0.0]))]
        zero=transient_field(P,zeros(N,1),[0.0]); n,dt=11,0.01
        force=transient_field(P,reshape(repeat([0.0,1.0],P.non),:,1),[0.0])
        exact,_=HHT(K,M,C,Global(force),zero,zero,n,dt;coordSys=cs)
        modified,_=HHT(K,M,C,zero,zero,zero,n,dt;coordSys=cs,
            modifyRHS=(b,s)->force)
        @test modified.a ≈ exact.a atol=1e-11
    end
    for options in ((reducedOrder=true,),(basis=:spectral,))
        with_test_gmsh() do
            line_mesh(n=2,order=3)
            P=Field([Material("body")];type=:ScalarField,dim=1,fieldName=:p,options...)
            M=∫(P⋅P); K=∫(Grad(P)⋅Grad(P))+M; N=ndofs(P)
            zero=transient_field(P,zeros(N,1),[0.0]); force=∫(P⋅1.0)
            support=[BoundaryCondition("left";field=P,p=(x,y,z,t)->0.1t)]
            exact,_=HHT(Symmetric(2K),Symmetric(M),force,zero,zero,13,0.01;support)
            modified,_=HHT(Symmetric(K),Symmetric(M),force,zero,zero,13,0.01;
                support,modifyA0=(A,s)->Symmetric(2K),solver=:cg,
                solveroptions=(reltol=1e-12,maxiter=200),copyOperators=false)
            @test modified.a ≈ exact.a atol=1e-10
        end
    end
end


@testset "modify operators: output validation" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        Q=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:q)
        N=ndofs(P); A=SystemMatrix(spdiagm(0=>ones(N)),P)
        z=transient_field(P,zeros(N,1),[0.0])
        @test_throws ErrorException HHT(A,A,z,z,z,3,0.1;
            modifyA0=(K,s)->SystemMatrix(copy(A.A),Q))
        @test_throws ErrorException HHT(A,A,z,z,z,3,0.1;
            modifyA0=(K,s)->SystemMatrix(spdiagm(0=>fill(NaN,N)),P))
        @test_throws DimensionMismatch HHT(A,A,z,z,z,3,0.1;
            modifyRHS=(b,s)->transient_field(P,zeros(N,2),[s.t,s.t+0.1]))
        @test_throws ArgumentError HHT(A,A,z,z,z,3,0.1;
            modifyA0=(K,s)->nothing,nonlinearoptions=(typo=1,))
        @test_throws ErrorException HHT(A,A,z,z,z,3,0.1;
            modifyA0=(K,s)->nothing,solver=:unknown)
    end
end

@testset "modify operators: private and borrowed matrix inputs" begin
    with_test_gmsh() do
        line_mesh(n=2, order=1)
        P = Field([Material("body")]; type=:ScalarField, dim=1, fieldName=:p)
        M = ∫(P⋅P)
        K = ∫(Grad(P)⋅Grad(P)) + M
        C = 0.2M
        z = transient_field(P, zeros(ndofs(P), 1), [0.0])
        load = ∫(P⋅1.0)
        saved = (copy(K.A), copy(M.A), copy(C.A), copy(load.a))
        n, dt = 6, 0.01
        support = [BoundaryCondition("left"; field=P, p=(x,y,z,t)->0.1t)]
        for run in (HHT, CDM, FDM), symmetric in (false, true)
            wrap(A) = symmetric ? Symmetric(A, :L) : A
            k, m, c = wrap(K), wrap(M), wrap(C)
            execute(; kwargs...) = run === FDM ?
                FDM(k, m, load, z, n, dt; support, kwargs...) :
                run(k, m, c, load, z, z, n, dt; support, kwargs...)
            result_fields(result) = result isa Tuple ? result : (result,)
            expected = result_fields(execute())
            for copyOperators in (true, false)
                seen = Set{Symbol}()
                function observe(key)
                    return (A, state) -> begin
                        accepted = getproperty(state.operators, key)
                        @test (A === accepted) == !copyOperators
                        @test (LowLevelFEM.parent_system(A).A ===
                               LowLevelFEM.parent_system(accepted).A) == !copyOperators
                        push!(seen, key)
                        # Exercise both unchanged-operator return conventions.
                        return key === :A0 ? A : nothing
                    end
                end
                function discard_rhs(rhs, state)
                    @test rhs.a !== load.a
                    rhs.a .= 99.0
                    return nothing
                end
                hooks = run === FDM ?
                    (; modifyA0=observe(:A0), modifyA1=observe(:A1)) :
                    (; modifyA0=observe(:A0), modifyA1=observe(:A1), modifyA2=observe(:A2))
                actual = result_fields(execute(; hooks..., modifyRHS=discard_rhs, copyOperators))
                for (a, b) in zip(actual, expected)
                    @test a.a ≈ b.a atol=1e-10
                end
                @test seen == (run === FDM ? Set((:A0, :A1)) : Set((:A0, :A1, :A2)))
            end
            # Default mode discards edits unless the private argument is returned.
            function discarded_edits(A, state)
                LowLevelFEM.parent_system(A).A .*= 2.0
                return nothing
            end
            actual = result_fields(execute(; modifyA0=discarded_edits))
            for (a, b) in zip(actual, expected)
                @test a.a ≈ b.a atol=1e-10
            end
            @test saved == (K.A, M.A, C.A, load.a)
            for invalid in (0, nothing, :off)
                @test_throws ArgumentError execute(; copyOperators=invalid)
                @test_throws ArgumentError execute(; copyOperators=invalid, modifyA0=(A,s)->nothing)
            end
        end
    end
end
