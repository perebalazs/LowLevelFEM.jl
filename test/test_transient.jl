using IterativeSolvers, Preconditioners, IncompleteLU

# Manufactured solutions exercise physical equilibrium, including moving BCs.
transient_field(P, A, t) = LowLevelFEM._transient_field(P, A, t)

@testset "Transient: damped oscillators and theta integration" begin
    with_test_gmsh() do
        line_mesh(n=1, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        N = ndofs(P)
        K = SystemMatrix(spdiagm(0=>fill(4.0,N)), P)
        M = SystemMatrix(spdiagm(0=>ones(N)), P)
        C = SystemMatrix(spdiagm(0=>fill(0.4,N)), P)
        zero = transient_field(P, zeros(N,1), [0.0])
        one = transient_field(P, ones(N,1), [0.0])
        n, dt = 401, 0.005
        t = collect(0:n-1)*dt
        exact = exp.(-0.2t).*(cos.(sqrt(3.96)*t) + 0.2/sqrt(3.96)*sin.(sqrt(3.96)*t))
        for run in (HHT, CDM)
            u,v = run(K,M,C,zero,one,zero,n,dt)
            @test maximum(abs.(u.a[1,:]-exact)) < (run === HHT ? 3e-5 : 0.002)
            @test u.t ≈ t
            @test u.kind === :field && v.kind === :field
            plain,_ = run(K,M,zero,one,zero,n,dt)
            nodamp,_ = run(K,M,0.0*C,zero,one,zero,n,dt)
            @test plain.a == nodamp.a
            kw,_ = run(K,M,zero,one,zero,n,dt; damping=C)
            @test kw.a == u.a
        end
        # Weighted force is essential for the exact linear response.
        ramp = transient_field(P, repeat(reshape(t,1,:),N,1), t)
        forcing = transient_field(P, repeat(reshape(4t .+ 0.4,1,:),N,1), t)
        u,v = HHT(Symmetric(K),Symmetric(M),Symmetric(C),forcing,zero,one,n,dt; α=-0.1)
        @test u.a ≈ ramp.a atol=2e-12
        @test v.a ≈ ones(N,n) atol=2e-12
        for theta in (0.0,0.5,1.0)
            q = transient_field(P, repeat(reshape(4t .+ 0.4,1,:),N,1), t)
            x = FDM(K,C,q,zero,n,dt; ϑ=theta)
            @test x.a ≈ ramp.a atol=2e-11
        end
        saved = copy(C.A),copy(forcing.a),copy(one.a)
        HHT(K,M,C,forcing,zero,one,n,dt)
        @test saved == (C.A,forcing.a,one.a)
        @test_throws ArgumentError HHT(K,M,zero,one,zero,1,dt)
        @test_throws ArgumentError CDM(K,M,zero,one,zero,n,NaN)
        @test_throws ArgumentError FDM(K,C,zero,one,n,dt;ϑ=1.1)
        wrong = transient_field(P,zeros(N,2),[0.0,1.0])
        @test_throws DimensionMismatch HHT(K,M,wrong,one,zero,n,dt)
        @test_throws DimensionMismatch HHT(K,M,zero,wrong,zero,n,dt)
    end
end

@testset "Transient: moving supports, inertia and damping coupling" begin
    with_test_gmsh() do
        line_mesh(n=1, order=1)
        P = Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        # Dense mass and damping cross terms catch omitted boundary derivatives.
        K = SystemMatrix(sparse([3.0 -1.0; -1.0 4.0]),P)
        M = SystemMatrix(sparse([2.0 0.3; 0.3 1.5]),P)
        C = SystemMatrix(sparse([0.4 0.1; 0.1 0.5]),P)
        n,dt = 101,0.01
        t = collect(0:n-1)*dt
        U = repeat(reshape(t.^2,1,:),2,1)
        V = repeat(reshape(2t,1,:),2,1)
        A = fill(2.0,2,n)
        bc_selector = [BoundaryCondition("left";field=P,p=0.0)]
        free = freeDoFs(P,bc_selector)
        U[free,:] .= 0.0
        V[free,:] .= 0.0
        A[free,:] .= 0.0
        initial = transient_field(P,zeros(2,1),[0.0])
        load = transient_field(P, M.A*A+C.A*V+K.A*U,t)
        support = [BoundaryCondition("left";field=P,p=(x,y,z,t)->t^2)]
        for run in (HHT,CDM)
            u,v = run(K,M,C,load,support,initial,initial,n,dt)
            @test u.a ≈ U atol=2e-11
            fixed = constrainedDoFs(P,support)
            @test v.a[fixed,:] ≈ V[fixed,:] atol=2e-10
            # A sampled BC includes support inertia without modifying initial fields.
            sampled_support=[BoundaryCondition("left";field=P,p=transient_field(P,U,t))]
            u2,_ = run(K,M,C,load,initial,initial,n,dt;support=sampled_support)
            @test u2.a ≈ U atol=2e-11
        end
        # First-order exact linear motion tests the capacity cross term.
        X = repeat(reshape(t,1,:),2,1)
        q = transient_field(P,C.A*ones(2,n)+K.A*X,t)
        bc = [BoundaryCondition("left";field=P,p=transient_field(P,X,t))]
        x = FDM(K,C,q,bc,initial,n,dt)
        @test x.a ≈ X atol=2e-12
        # Spatial-only callbacks and a shifted physical time grid.
        bcall = [BoundaryCondition("body";field=P,p=(x,y,z,t)->t^2)]
        u,_ = HHT(K,M,load,initial,initial,n,dt;support=bcall,t0=2.0)
        @test u.a ≈ repeat(reshape((t .+ 2).^2,1,:),2,1)
        @test u.t ≈ t .+ 2
        bcall = [BoundaryCondition("body";field=P,p=(x,y,z)->3.0)]
        x = FDM(K,C,initial,initial,2,dt;support=bcall)
        @test x.a == fill(3.0,2,2)
    end
end

@testset "Transient: accuracy analysis with physical damping" begin
    dt = 0.1
    _,undamped = HHTaccuracyAnalysis(1.0,3.0,dt,:SR;n=5)
    _,damped = HHTaccuracyAnalysis(1.0,3.0,dt,:SR;n=5,ξ=0.05)
    @test all(damped .< undamped)
    _,ratios = HHTaccuracyAnalysis(1.0,3.0,dt,:PDR;n=5,ξ=0.02,massDamping=0.1,stiffnessDamping=0.04)
    omega = collect(range(1.0,3.0;length=5))
    @test ratios ≈ 0.02 .+ 0.1 ./ (2omega) .+ 0.02omega
    _,radius = FDMaccuracyAnalysis(2.0,2.0,dt;n=2,ϑ=1.0)
    @test radius ≈ fill(1/1.2,2)
    @test length(first(CDMaccuracyAnalysis(1.0,3.0,dt,:SR))) == 100
end

@testset "Transient: multifield dynamics, MPC and block ordering" begin
    with_test_gmsh() do
        line_mesh(n=2,order=1)
        mat = Material("body")
        P = Field([mat],type=:ScalarField,dim=1,fieldName=:p)
        S = Field([mat],type=:ScalarField,dim=1,fieldName=:s)
        N = ndofs(P)
        Ip = SystemMatrix(spdiagm(0=>ones(N)),P)
        Is = SystemMatrix(spdiagm(0=>ones(N)),S)
        Z = SystemMatrix(spzeros(N,N),S,P)
        K = SystemMatrix([2Ip Z; Z' 3Is])
        M = SystemMatrix([Ip Z; Z' Is])
        C = 0.2M
        n,dt=51,0.01
        t=collect(0:n-1)*dt
        up = repeat(reshape(t,1,:),N,1)
        us = 2up
        zeroP,zeroS = transient_field(P,zeros(N,1),[0.0]), transient_field(S,zeros(N,1),[0.0])
        vP,vS = transient_field(P,ones(N,1),[0.0]), transient_field(S,fill(2.0,N,1),[0.0])
        qp = transient_field(P,2up .+ 0.2,t)
        qs = transient_field(S,3us .+ 0.4,t)
        F = SystemVector([qp,qs])
        U0,V0 = SystemVector([zeroP,zeroS]),SystemVector([vP,vS])
        mpc = [MPC(master="left",slave="right",field=P,p=true)]
        for run in (HHT,CDM)
            (p,s),(vp,vs) = run(K,M,C,F,U0,V0,n,dt;mpc)
            @test p.a ≈ up atol=2e-12
            @test s.a ≈ us atol=2e-12
            @test vp.a ≈ ones(N,n) atol=2e-11
            @test vs.a ≈ fill(2.0,N,n) atol=2e-11
            @test p.t == s.t == t
            # Field tuples are accepted in block order as well.
            (p2,s2),_ = run(K,M,C,(qp,qs),(zeroP,zeroS),(vP,vS),n,dt;mpc)
            @test p2.a ≈ p.a atol=2e-12
            @test s2.a ≈ s.a atol=2e-12
        end
        p,s = FDM(K,C,F,U0,n,dt;mpc)
        @test p.a ≈ up atol=2e-12
        @test s.a ≈ us atol=2e-12
        # Moving slave constraints must transfer to the master for every sample.
        support=[BoundaryCondition("right";field=P,p=(x,y,z,t)->t)]
        (p,s),_ = HHT(K,M,C,F,U0,V0,n,dt;mpc,support)
        @test p.a ≈ up atol=2e-12
        @test s.a ≈ us atol=2e-12
        supportbad=[BoundaryCondition("left";field=P,p=0.0),
                    BoundaryCondition("right";field=P,p=(x,y,z,t)->t)]
        @test_throws ErrorException HHT(K,M,C,F,U0,V0,n,dt;mpc,support=supportbad)
        wrong=SystemVector([qs,qp])
        @test_throws ErrorException FDM(K,C,wrong,U0,n,dt)
        @test_throws ErrorException FDM(K,C,F,U0,n,dt;support=[BoundaryCondition("left";p=0.0)])
    end
end

@testset "Transient: local coordinates and sampled vector supports" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1,order=1)
        P=Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        N=ndofs(P)
        K=SystemMatrix(spdiagm(0=>fill(2.0,N)),P)
        M=SystemMatrix(spdiagm(0=>ones(N)),P)
        C=0.3M
        e1=VectorField(P,"body",[0.0,1.0,0.0])
        cs=CoordinateSystem(e1)
        Q=LowLevelFEM._build_coordinate_transformation(P,[cs]).T
        n,dt=21,0.01
        t=collect(0:n-1)*dt
        speed=repeat([1.0,2.0],P.non)
        U=speed*t'
        V=repeat(reshape(speed,:,1),1,n)
        initial=transient_field(P,zeros(N,1),[0.0])
        velocity=transient_field(P,reshape(speed,:,1),[0.0])
        f=transient_field(P,2U+0.3V,t)
        fg=transient_field(P,Q*f.a,t)
        vg=transient_field(P,Q*velocity.a,[0.0])
        # Only the local x component is constrained; y remains dynamic.
        support=[BoundaryCondition("left";field=P,ux=ScalarField([],repeat(reshape(t,1,:),P.non,1),t,[],n,:scalar,P))]
        d=transient_field(P,U,t)
        for run in (HHT,CDM)
            u,v=run(K,M,C,f,initial,velocity,n,dt;coordSys=[cs],support)
            @test u.a ≈ Q*U atol=2e-12
            @test v.a ≈ Q*V atol=2e-11
            globalu,_=run(K,M,C,Global(fg),Global(initial),Global(vg),n,dt;
                           coordSys=[cs],support)
            @test globalu.a ≈ u.a atol=2e-12
        end
        x=FDM(K,C,Global(fg),initial,n,dt;coordSys=[cs],support)
        @test x.a ≈ Q*U atol=2e-12
        # Local and Cartesian load terms can be combined columnwise.
        constant=transient_field(P,fill(0.1,N,n),t)
        remainder=transient_field(P,Q*(f.a .- constant.a),t)
        mixed,_=HHT(K,M,C,constant+Global(remainder),initial,velocity,n,dt;coordSys=[cs])
        @test mixed.a ≈ Q*U atol=2e-12
    end
end

@testset "Transient: reducedOrder and SEM projection" begin
    for options in ((reducedOrder=true,), (basis=:spectral,))
        with_test_gmsh() do
            line_mesh(n=2,order=3)
            P=Field([Material("body")];type=:ScalarField,dim=1,fieldName=:p,options...)
            K=∫(Grad(P)⋅Grad(P))+∫(P⋅P)
            M=∫(P⋅P)
            C=0.1M
            N=ndofs(P)
            n,dt=21,0.002
            t=collect(0:n-1)*dt
            U=repeat(reshape(t,1,:),N,1)
            V=ones(N,n)
            f=transient_field(P,K.A*U+C.A*V,t)
            zero=transient_field(P,zeros(N,1),[0.0])
            speed=transient_field(P,ones(N,1),[0.0])
            support=[BoundaryCondition("left";field=P,p=(x,y,z,t)->t)]
            for run in (HHT,CDM)
                u,v=run(K,M,C,f,zero,speed,n,dt;support)
                @test u.a ≈ U atol=2e-11
                @test v.a ≈ V atol=2e-10
            end
            x=FDM(K,C,f,zero,n,dt;support)
            @test x.a ≈ U atol=2e-11
            iterative,_=HHT(Symmetric(K),Symmetric(M),Symmetric(C),f,zero,speed,n,dt;
                            support,solver=:cg,solveroptions=(reltol=1e-12,maxiter=200))
            @test iterative.a ≈ U atol=2e-10
            iterative_x=FDM(K,C,f,zero,n,dt;support,solver=:gmres,
                            solveroptions=(reltol=1e-12,maxiter=200))
            @test iterative_x.a ≈ U atol=2e-10
            # An independent reduced-system reference verifies nonconstant modes.
            T,R=LowLevelFEM._approximation_matrices(P)
            kr,mr,cr=T'*K.A*T,T'*M.A*T,T'*C.A*T
            r0=collect(1.0:size(T,2))/size(T,2)
            u0=transient_field(P,reshape(T*r0,:,1),[0.0])
            f0=transient_field(P,zeros(N,1),[0.0])
            u,v=HHT(K,M,C,f0,u0,zero,n,dt)
            a=mr\(-kr*r0)
            r=copy(r0); vel=zeros(length(r))
            effective=mr+dt/2*cr+dt^2/4*kr
            for j in 2:n
                pred=r+dt*vel+dt^2/4*a
                vpred=vel+dt/2*a
                a=effective\(-kr*pred-cr*vpred)
                r=pred+dt^2/4*a
                vel=vpred+dt/2*a
                @test u.a[:,j] ≈ T*r atol=2e-11
                @test v.a[:,j] ≈ T*vel atol=2e-10
            end
        end
    end
end

@testset "Transient: reduced scalar + rotated vector blocks with MPC" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1,order=2)
        mat=Material("body")
        P=Field([mat],type=:ScalarField,dim=2,fieldName=:p,reducedOrder=true)
        V=Field([mat],type=:VectorField,dim=2,fieldName=:v)
        Mp,Mv=∫(P⋅P),∫(V⋅V)
        Z=SystemMatrix(spzeros(ndofs(P),ndofs(V)),V,P)
        M=SystemMatrix([Mp Z; Z' Mv])
        K=SystemMatrix([2Mp Z; Z' 3Mv])
        C=0.2M
        cs=[CoordinateSystem(VectorField(V,"body",[0.0,1.0,0.0]))]
        Q=LowLevelFEM._build_multifield_coordinate_transformation(K,cs)
        n,dt=11,0.01
        t=collect(0:n-1)*dt
        speed=vcat(ones(ndofs(P)),repeat([1.0,2.0],V.non))
        U=speed*t'
        globalF=K.A*(Q*U)+C.A*(Q*repeat(reshape(speed,:,1),1,n))
        block(A)=SystemVector(Matrix(A),nothing,K.problems,K.offsets)
        F=Global(block(globalF))
        U0=block(zeros(size(M.A,1),1))
        V0=block(reshape(speed,:,1))
        ties=[MPC(master="leftbottom",slave="rightbottom",field=V,vx=true,vy=true)]
        bc=[BoundaryCondition("left";field=P,p=(x,y,z,t)->t),
            BoundaryCondition("top";field=V,vx=(x,y,z,t)->t)]
        for run in (HHT,CDM)
            (p,v),_=run(K,M,C,F,U0,V0,n,dt;support=bc,mpc=ties,coordSys=cs)
            @test vcat(p.a,v.a) ≈ Q*U atol=2e-11
        end
        p,v=FDM(K,C,F,U0,n,dt;support=bc,mpc=ties,coordSys=cs)
        @test vcat(p.a,v.a) ≈ Q*U atol=2e-11
        (p,v),_=HHT(Symmetric(K),Symmetric(M),Symmetric(C),F,U0,V0,n,dt;
                    support=bc,mpc=ties,coordSys=cs,solver=:cg,
                    solveroptions=(reltol=1e-12,maxiter=200))
        @test vcat(p.a,v.a) ≈ Q*U atol=2e-10
        p,v=FDM(K,C,F,U0,n,dt;support=bc,mpc=ties,coordSys=cs,solver=:gmres,
                 solveroptions=(reltol=1e-12,maxiter=200))
        @test vcat(p.a,v.a) ≈ Q*U atol=2e-10
    end
end

@testset "Transient: sampled scalar BCs, tensor fields and validation" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        I2=SystemMatrix(spdiagm(0=>ones(ndofs(P))),P)
        zero=transient_field(P,zeros(ndofs(P),1),[0.0])
        n,dt=5,0.05
        t=collect(0:n-1)*dt
        sampled=transient_field(P,repeat(reshape(t,1,:),ndofs(P),1),t)
        bc=[BoundaryCondition("body";field=P,p=sampled)]
        u,v=HHT(I2,I2,zero,bc,zero,zero,n,dt)
        @test u.a == sampled.a
        @test v.a ≈ ones(ndofs(P),n) atol=1e-12
        u,v=CDM(I2,I2,zero,bc,zero,zero,n,dt)
        @test u.a == sampled.a
        @test v.a ≈ ones(ndofs(P),n) atol=1e-12
        @test_throws DimensionMismatch HHT(I2,I2,zero,sampled,zero,n,dt)
        @test_throws ErrorException HHT(I2,I2,zero,zero,zero,n,dt;
                    support=[BoundaryCondition("left";field=P,p=NaN)])
        # Legacy names and zero-damping positional signatures are retained.
        @test LowLevelFEM.FDM_old(I2,I2,zero,BoundaryCondition[],zero,n,dt).a ==
              FDM(I2,I2,zero,zero,n,dt).a
    end
    with_test_gmsh() do
        structured_rect_mesh(n=1,order=1)
        P=Field([Material("body")],type=:TensorField,dim=2,fieldName=:s)
        N=ndofs(P)
        M=SystemMatrix(spdiagm(0=>ones(N)),P)
        K=2M
        C=0.2M
        n,dt=5,0.01
        t=collect(0:n-1)*dt
        zero=transient_field(P,zeros(N,1),[0.0])
        speed=transient_field(P,ones(N,1),[0.0])
        expected=repeat(reshape(t,1,:),N,1)
        f=transient_field(P,2expected .+ 0.2,t)
        for run in (HHT,CDM)
            u,v=run(K,M,C,f,zero,speed,n,dt)
            @test u isa TensorField && v isa TensorField
            @test u.a ≈ expected atol=1e-12
        end
        @test FDM(K,C,f,zero,n,dt).a ≈ expected atol=1e-12
    end
end

@testset "HHT: weighted equilibrium with nonproportional damping" begin
    with_test_gmsh() do
        line_mesh(n=2,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        K=SystemMatrix(sparse([3.0 -1.0 0; -1 4 -1; 0 -1 2]),P)
        M=SystemMatrix(sparse([2.0 0.2 0; 0.2 1.5 0.1; 0 0.1 1]),P)
        C=SystemMatrix(sparse([0.4 0.1 0; 0.1 0.8 0.2; 0 0.2 0.3]),P)
        n,dt=31,0.01
        t=collect(0:n-1)*dt
        force=repeat(reshape(sin.(3t),1,:),ndofs(P),1)
        initial=transient_field(P,zeros(ndofs(P),1),[0.0])
        D=repeat(reshape(sin.(t),1,:),ndofs(P),1)
        support=[BoundaryCondition("left";field=P,p=transient_field(P,D,t))]
        free,fixed=freeDoFs(P,support),constrainedDoFs(P,support)
        alpha,gamma,beta=-0.1,0.6,0.3025
        u,v=HHT(K,M,C,transient_field(P,force,t),initial,initial,n,dt;
                 α=alpha,γ=gamma,β=beta,support)
        a=zeros(ndofs(P),n)
        a[fixed,1].=(2D[fixed,1]-5D[fixed,2]+4D[fixed,3]-D[fixed,4])/dt^2
        a[free,1].=M.A[free,free]\(force[free,1]-K.A[free,:]*u.a[:,1]-
                                   C.A[free,:]*v.a[:,1]-M.A[free,fixed]*a[fixed,1])
        for j in 2:n
            a[:,j].=(v.a[:,j]-v.a[:,j-1]-dt*(1-gamma)*a[:,j-1])/(gamma*dt)
            equilibrium=M.A*a[:,j]+(1+alpha)*(K.A*u.a[:,j]+C.A*v.a[:,j]-force[:,j])-
                        alpha*(K.A*u.a[:,j-1]+C.A*v.a[:,j-1]-force[:,j-1])
            @test norm(equilibrium[free]) < 2e-11
            @test u.a[free,j] ≈ u.a[free,j-1]+dt*v.a[free,j-1]+
                  dt^2*((0.5-beta)*a[free,j-1]+beta*a[free,j]) atol=2e-12
        end
    end
end

@testset "HHT: lumped mass and coupled stiffness keep dense states" begin
    with_test_gmsh() do
        structured_rect_mesh(n=2)
        mat = Material("body", E=2e5, ν=0.3, ρ=7.86e-9)
        P = Problem([mat], type=:PlaneStress)
        K = stiffnessMatrix(P)
        M = massMatrix(P; lumped=true)
        support = [displacementConstraint("left"; ux=0.0, uy=0.0)]
        u0 = initialDisplacement(P, "body"; ux=0.0, uy=0.0)
        v0 = initialVelocity(P, "body"; vx=0.0, vy=0.0)
        n, dt = 21, 1e-7
        applied(x,y,z,t) = t < 4dt ? t/(4dt) : t < 8dt ? (8dt-t)/(4dt) : 0.0
        fx = ScalarField(P, "right", applied; steps=n, tmax=(n-1)*dt)
        f = loadVector(P, [load("right"; fx=fx)])
        free = freeDoFs(P, support)
        fixed = constrainedDoFs(P, support)
        @test isdiag(M.A)
        @test !isdiag(K.A[free, free])
        for C in (nothing, 0.01M)
            u,v = HHT(K,M,f,u0,v0,n,dt; support, damping=C)
            @test all(isfinite, u.a) && all(isfinite, v.a)
            @test iszero(u.a[fixed,:]) && iszero(v.a[fixed,:])
            @test !iszero(u.a[free,:])
            # Independent dense average-acceleration Newmark recurrence.
            k, m = Matrix(K.A[free,free]), Matrix(M.A[free,free])
            c = C === nothing ? zeros(size(m)) : Matrix(C.A[free,free])
            x, vel = zeros(length(free)), zeros(length(free))
            a = m \ f.a[free,1]
            effective = m + dt/2*c + dt^2/4*k
            for j in 2:n
                xp = x + dt*vel + dt^2/4*a
                vp = vel + dt/2*a
                a = effective \ (f.a[free,j] - k*xp - c*vp)
                x = xp + dt^2/4*a
                vel = vp + dt/2*a
                @test u.a[free,j] ≈ x rtol=1e-10 atol=1e-14
                @test v.a[free,j] ≈ vel rtol=1e-10 atol=1e-12
            end
        end
    end
end

@testset "HHT/FDM: reusable direct and iterative solvers" begin
    with_test_gmsh() do
        line_mesh(n=8, order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        M=∫(P⋅P)
        K=∫(Grad(P)⋅Grad(P))+M
        C=0.3M
        n,dt=17,0.01
        t=collect(0:n-1)*dt
        N=ndofs(P)
        spatial=collect(1.0:N)/N
        f=transient_field(P,spatial*(1 .+ sin.(5t))',t)
        zero=transient_field(P,zeros(N,1),[0.0])
        support=[BoundaryCondition("left";field=P,p=(x,y,z,t)->0.01sin(2t))]
        reference,velocity=HHT(K,M,C,f,zero,zero,n,dt;support,α=-0.1,δ=0.1,solver=:lu)
        diffusion=FDM(K,M,f,zero,n,dt;support,ϑ=0.7,solver=:lu)
        wrapped=(Symmetric(K),Symmetric(M),Symmetric(C))
        for solver in (:auto,:backslash,:lu,:lu_no_ordering,:cholesky,:ldlt,:qr,:cg,:gmres)
            options=solver in (:cg,:gmres) ? (reltol=1e-12,maxiter=200) : (;)
            u,v=HHT(wrapped...,f,zero,zero,n,dt;support,α=-0.1,δ=0.1,
                     solver,solveroptions=options)
            x=FDM(wrapped[1],wrapped[2],f,zero,n,dt;support,ϑ=0.7,
                  solver,solveroptions=options)
            @test relerr(u.a,reference.a) < 1e-9
            @test relerr(v.a,velocity.a) < 1e-8
            @test relerr(x.a,diffusion.a) < 1e-9
        end
        free=freeDoFs(P,support)
        pc=Diagonal(Vector(diag((M.A+dt^2/4*K.A)[free,free])))
        for (solver,preconditioner) in ((:cg,:none),(:cg,:ich),(:gmres,:none),
                                       (:gmres,:ilu),(:cg,pc),(:gmres,pc))
            u,v=HHT(wrapped...,f,zero,zero,n,dt;support,α=-0.1,δ=0.1,
                     solver,preconditioner,solveroptions=(reltol=1e-12,maxiter=200))
            x=FDM(wrapped[1],wrapped[2],f,zero,n,dt;support,ϑ=0.7,
                  solver,preconditioner,solveroptions=(reltol=1e-12,maxiter=200))
            @test relerr(u.a,reference.a) < 1e-9
            @test relerr(v.a,velocity.a) < 1e-8
            @test relerr(x.a,diffusion.a) < 1e-9
        end
        # Native preconditioner options follow solveField without rebuilding per step.
        u,_=HHT(wrapped...,f,zero,zero,n,dt;support,α=-0.1,δ=0.1,solver=:cg,
                 preconditioneroptions=(memory=3,),solveroptions=(reltol=1e-12,))
        x=FDM(K,M,f,zero,n,dt;support,ϑ=0.7,solver=:gmres,
              preconditioneroptions=(τ=1e-3,),solveroptions=(reltol=1e-12,))
        @test relerr(u.a,reference.a) < 1e-9
        @test relerr(x.a,diffusion.a) < 1e-9
        for opts in ((solver=:unknown,), (solver=:lu,solveroptions=(reltol=1e-6,)),
                     (solver=:lu,preconditioner=:ilu),
                     (solver=:lu,preconditioneroptions=(τ=0.01,)),
                     (solver=:gmres,solveroptions=(log=true,)),
                     (solver=:gmres,solveroptions=[1,2]),
                     (solver=:gmres,preconditioner=:bad),
                     (solver=:gmres,preconditioner=:none,preconditioneroptions=(τ=0.01,)),
                     (solver=:cg,preconditioner=pc,preconditioneroptions=(memory=3,)))
            @test_throws ErrorException HHT(K,M,C,f,zero,zero,n,dt;opts...)
            @test_throws ErrorException FDM(K,M,f,zero,n,dt;opts...)
        end
        # A nonsymmetric damping matrix must not inherit K/M symmetry promises.
        nonsymmetric=SystemMatrix(C.A+spdiagm(1=>fill(0.01,N-1)),P)
        u,_=HHT(wrapped[1],wrapped[2],nonsymmetric,f,zero,zero,n,dt;support,solver=:gmres,
                 solveroptions=(reltol=1e-12,))
        exact,_=HHT(K,M,nonsymmetric,f,zero,zero,n,dt;support,solver=:lu)
        @test relerr(u.a,exact.a) < 1e-9
        @test_throws ErrorException HHT(wrapped[1],wrapped[2],nonsymmetric,f,zero,zero,n,dt;
                                        support,solver=:cg)
    end
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        M=SystemMatrix(spdiagm(0=>ones(ndofs(P))),P)
        K=-2M
        zero=transient_field(P,zeros(ndofs(P),1),[0.0])
        one=transient_field(P,ones(ndofs(P),1),[0.0])
        @test_logs (:warn,r"Cholesky factorization failed") begin
            x=FDM(Symmetric(K),Symmetric(M),zero,one,2,1.0;ϑ=1.0)
            @test x.a[:,2] ≈ -ones(ndofs(P))
        end
        @test_throws PosDefException FDM(Symmetric(K),Symmetric(M),zero,one,2,1.0;
                                         ϑ=1.0,solver=:cholesky)
    end
end


@testset "Transient: support is the only prescribed-history input" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P=Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        N=ndofs(P); K=SystemMatrix(spdiagm(0=>fill(2.0,N)),P)
        M=SystemMatrix(spdiagm(0=>ones(N)),P); C=0.2M
        n,dt=7,0.01; t=collect(0:n-1)*dt
        zero=transient_field(P,zeros(N,1),[0.0])
        history=transient_field(P,repeat(reshape(t,1,:),N,1),t)
        support=[BoundaryCondition("left";field=P,p=history)]
        fixed=constrainedDoFs(P,support)
        force=transient_field(P,repeat(reshape(2t .+ 0.2,1,:),N,1),t)
        wrong_velocity=transient_field(P,fill(99.0,N,1),[0.0])
        for run in (HHT,CDM)
            u,v=run(K,M,C,force,zero,wrong_velocity,n,dt;support)
            @test u.a[fixed,:] ≈ reshape(t,1,:)
            @test v.a[fixed,:] ≈ ones(1,n) atol=1e-11
            @test_throws DimensionMismatch run(K,M,zero,history,zero,n,dt;support)
        end
        x=FDM(K,C,force,zero,n,dt;support)
        @test x.a[fixed,:] ≈ reshape(t,1,:)
        @test_throws DimensionMismatch FDM(K,C,force,history,n,dt;support)
        @test iszero(zero.a)
        @test history.a == repeat(reshape(t,1,:),N,1)
        short=transient_field(P,zeros(N,n-1),t[1:end-1])
        @test_throws ErrorException HHT(K,M,zero,zero,zero,n,dt;
            support=[BoundaryCondition("left";field=P,p=short)])
    end
end
