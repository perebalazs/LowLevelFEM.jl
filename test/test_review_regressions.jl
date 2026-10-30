@testset "Review: HHT moving-support convergence with consistent mass" begin
    with_test_gmsh() do
        line_mesh(n=1, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        K = SystemMatrix(sparse([3.0 -1.0; -1.0 4.0]), P)
        M = SystemMatrix(sparse([2.0 0.6; 0.6 1.5]), P)
        C = SystemMatrix(sparse([0.4 0.1; 0.1 0.5]), P)
        z = transient_field(P, zeros(2,1), [0.0])
        selector = [BoundaryCondition("left"; field=P, p=0.0)]
        f, c = freeDoFs(P, selector), constrainedDoFs(P, selector)
        for alpha in (0.0, -0.1)
            errors = Float64[]
            for steps in (201, 401, 801)
                dt = 1/(steps-1)
                t = collect(0:steps-1)*dt
                D = zeros(2,steps); D[c,:] .= sin.(2π*t)'
                V = zeros(2,steps); V[c,:] .= (2π*cos.(2π*t))'
                A = -(2π)^2 .* D
                load = transient_field(P, M.A*A + C.A*V + K.A*D, t)
                support = [BoundaryCondition("left"; field=P, p=(x,y,z,t)->sin(2π*t))]
                u,v = HHT(K,M,C,load,z,z,steps,dt; support, α=alpha)
                acceleration_error = Ref(0.0)
                function observe(rhs, state)
                    actual = state.acceleration[1].a[c[1],1]
                    exact = -(2π)^2*sin(2π*state.t)
                    acceleration_error[] = max(acceleration_error[],abs(actual-exact))
                    return nothing
                end
                um,vm = HHT(K,M,C,load,z,z,steps,dt; support, α=alpha,
                            modifyRHS=observe,nonlinearoptions=(reltol=1e-12,abstol=1e-13))
                @test acceleration_error[] < 0.01
                @test um.a ≈ u.a atol=2e-10
                @test vm.a ≈ v.a atol=2e-9
                @test maximum(abs.(v.a[c,:]-V[c,:])) < 0.01
                push!(errors, maximum(abs, u.a[f,:]))
            end
            @test errors[1]/errors[2] > 3.5
            @test errors[2]/errors[3] > 3.5
            @test last(errors) < 3e-5
        end
        initial = transient_field(P,ones(2,1),[0.0])
        auto,_ = HHT(K,M,C,z,initial,z,13,0.01; α=-0.1)
        explicit,_ = HHT(K,M,C,z,initial,z,13,0.01; α=-0.1,γ=0.6,β=0.3025)
        historical,_ = HHT(K,M,C,z,initial,z,13,0.01; α=-0.1,δ=0.0)
        @test auto.a ≈ explicit.a
        @test !isapprox(auto.a,historical.a;rtol=1e-9)
        @test HHTaccuracyAnalysis(1.,3.,0.1,:SR;α=-0.1) ==
              HHTaccuracyAnalysis(1.,3.,0.1,:SR;α=-0.1,δ=0.1)
    end
end

@testset "Review: compact boundary storage and callback signatures" begin
    with_test_gmsh() do
        line_mesh(n=100, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        N = ndofs(P)
        I0 = SystemMatrix(spdiagm(0=>ones(N)),P)
        z = transient_field(P,zeros(N,1),[0.0])
        support = [BoundaryCondition("left"; field=P, p=(x,y,z,t)->sin(t))]
        prepare() = LowLevelFEM.prepare_transient_system(I0,I0,nothing,z,z,z,4001,0.001;support)
        p = prepare()
        @test size(p.D) == (1,4001)
        @test p.D[1,:] ≈ sin.(p.t)
        # A single old full history here required >3 MB, before its copies.
        @test (@allocated prepare()) < 3_000_000
        varargs(x,y,z,args...) = isempty(args) ? 3.0 : 99.0
        optional(x,y,z,t=7.0) = t
        for (callback, expected) in ((varargs,3.0),(optional,7.0))
            x = FDM(I0,I0,z,z,4,0.1; support=[BoundaryCondition("left";field=P,p=callback)])
            @test all(==(expected), x.a[constrainedDoFs(P,support),:])
        end
    end
end

@testset "Review: authoritative symmetric triangles before projection" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1,order=1)
        P = Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        N = ndofs(P)
        K = SystemMatrix(spdiagm(-1=>fill(-0.2,N-1),0=>collect(2.0:N+1),1=>fill(-0.2,N-1)),P)
        M = SystemMatrix(spdiagm(-1=>fill(0.03,N-1),0=>ones(N),1=>fill(0.03,N-1)),P)
        C = 0.1M
        z = transient_field(P,zeros(N,1),[0.0])
        load = transient_field(P,ones(N,1),[0.0])
        coordSys = [CoordinateSystem(VectorField(P,"body",[0.0,1.0,0.0]))]
        mpc = [MPC(master="left",slave="right",field=P,uy=true)]
        support = [BoundaryCondition("left";field=P,ux=0.0)]
        opts = (;support,mpc,coordSys)
        full = solveField(K,load;opts...)
        λ,ϕ = solveEigen(K,M;n=2,solver=:direct,normalize=:mass,output=:raw,opts...)
        uh,vh = HHT(K,M,C,load,z,z,13,0.01;opts...)
        uc,vc = CDM(K,M,C,load,z,z,13,0.01;opts...)
        xf = FDM(K,M,load,z,13,0.01;opts...)
        for (triangle,uplo) in ((triu,:U),(tril,:L))
            wrap(A) = Symmetric(SystemMatrix(triangle(A.A),P),uplo)
            k,m,c = wrap(K),wrap(M),wrap(C)
            @test solveField(k,load;opts...).a ≈ full.a
            vals,modes = solveEigen(k,m;n=2,solver=:direct,normalize=:mass,output=:raw,opts...)
            @test vals ≈ λ
            @test modes ≈ ϕ atol=1e-10
            @test modes'*M.A*modes ≈ Matrix{Float64}(I,2,2) atol=1e-11
            for (run,expected) in ((HHT,uh),(CDM,uc))
                u,_ = run(k,m,c,load,z,z,13,0.01;opts...)
                um,_ = run(k,m,c,load,z,z,13,0.01;opts...,modifyA0=(A,s)->A)
                @test u.a ≈ expected.a atol=1e-12
                @test um.a ≈ expected.a atol=1e-11
                changed,_ = run(k,m,c,load,z,z,13,0.01;opts...,
                    modifyA0=(A,s)->wrap(2K),modifyA2=(A,s)->wrap(3M),modifyA1=(A,s)->wrap(2C))
                reference,_ = run(2K,3M,2C,load,z,z,13,0.01;opts...)
                @test changed.a ≈ reference.a atol=1e-11
            end
            @test FDM(k,m,load,z,13,0.01;opts...).a ≈ xf.a
            updated = FDM(k,m,load,z,13,0.01;opts...,
                modifyA0=(A,s)->wrap(2K),modifyA1=(A,s)->wrap(3M))
            @test updated.a ≈ FDM(2K,3M,load,z,13,0.01;opts...).a atol=1e-11
        end
    end
    @test_throws SingularException LowLevelFEM.prepare_linear_solver(spdiagm(0=>[1.0,0.0]))
end

@testset "Review: centered CDM damping and velocity output" begin
    with_test_gmsh() do
        line_mesh(n=1,order=1)
        P = Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        N = ndofs(P)
        I0 = SystemMatrix(spdiagm(0=>ones(N)),P)
        z = transient_field(P,zeros(N,1),[0.0])
        u0 = transient_field(P,ones(N,1),[0.0])
        errors = Float64[]
        for dt in (0.02,0.01,0.005)
            n = round(Int,2/dt)+1; t = collect(0:n-1)*dt
            ω = sqrt(3.96)
            exact = exp.(-0.2t).*(cos.(ω*t)+0.2/ω*sin.(ω*t))
            speed = -4/ω*exp.(-0.2t).*sin.(ω*t)
            options = (;dampingScheme=:central,velocityScheme=:central)
            u,v = CDM(4I0,I0,0.4I0,z,u0,z,n,dt;options...)
            um,vm = CDM(4I0,I0,0.4I0,z,u0,z,n,dt;options...,modifyA0=(A,s)->nothing)
            @test um.a ≈ u.a atol=1e-11
            @test vm.a ≈ v.a atol=1e-10
            push!(errors,max(maximum(abs.(u.a[1,:]-exact)),maximum(abs.(v.a[1,:]-speed))))
        end
        @test errors[1]/errors[2] > 3.7
        @test errors[2]/errors[3] > 3.7
        @test_throws ArgumentError CDM(I0,I0,z,z,z,3,0.1;dampingScheme=:wrong)
        Ω,ξ = 0.2,0.1
        _,ρ = CDMaccuracyAnalysis(1.,1.,Ω,:SR;n=2,ξ,dampingScheme=:central)
        @test ρ[1] ≈ maximum(abs,eigvals([(2-Ω^2)/(1+ξ*Ω) (ξ*Ω-1)/(1+ξ*Ω);1 0]))
    end
end

@testset "Review: eigen selection, phase, and dense memory guard" begin
    with_test_gmsh() do
        line_mesh(n=39,order=1)
        P = Field([Material("body")],type=:ScalarField,dim=1,fieldName=:p)
        N = ndofs(P)
        values = vcat(0.4999999,0.50000001,collect(2.0:N-2),1e8)
        K = SystemMatrix(spdiagm(0=>values),P)
        M = SystemMatrix(spdiagm(0=>ones(N)),P)
        expected,_ = solveEigen(K,M;n=1,shift=0.5,solver=:direct,output=:raw)
        actual,modes = solveEigen(K,M;n=1,shift=0.5,solver=:arpack,output=:raw,
                                 solveroptions=(tol=1e-12,ncv=20,v0=ones(N)))
        @test actual ≈ expected atol=1e-10
        pivot = modes[argmax(abs.(modes[:,1])),1]
        @test real(pivot)>0 && abs(imag(pivot))<1e-12
        @test_throws ArgumentError solveEigen(K,M;n=2,solver=:direct,dense_limit=100)
        @test_throws ArgumentError solveEigen(K,M;n=2,solver=:direct,dense_limit=0)
        @test length(first(solveEigen(K,M;n=2,solver=:direct,dense_limit=nothing,output=:raw)))==2
    end
end

# Controlled failures exercise exception propagation without allocating memory
# or interrupting the test process itself.
struct _ReviewFailingMetric <: AbstractMatrix{Float64}
    failure::Exception
end
Base.size(::_ReviewFailingMetric) = (2,2)
Base.getindex(::_ReviewFailingMetric,i::Int,j::Int) = Float64(i==j)
LinearAlgebra.cholesky(B::Symmetric{Float64,_ReviewFailingMetric}) = throw(parent(B).failure)

@testset "Review: eigen classification propagates unexpected failures" begin
    for failure in (InterruptException(),OutOfMemoryError(),ArgumentError("injected"))
        prepared = (;A=spdiagm(0=>[1.0,2.0]),B=_ReviewFailingMetric(failure),
                     symmetricA=true,symmetricB=true,uploA=:U,uploB=:U)
        @test_throws typeof(failure) LowLevelFEM.solve_eigen_reduced(prepared;n=1,solver=:direct)
    end
end
