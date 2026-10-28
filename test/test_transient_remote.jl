# A remote cross-section is described by translations u and rotations φ.
# Only the displacement material carries physical mass; rotational inertia
# and translation/rotation coupling must follow from the same rigid map.
function transient_remote_fixture(dim)
    if dim == 2
        structured_rect_mesh(n=2,rightRemote=true)
    else
        structured_box_mesh(n=2,rightRemote=true)
    end
    mat=Material("body",E=100.0,ν=0.25,ρ=1.0)
    materials=[mat,Material("rightremote")]
    U=Field(materials,type=:VectorField,dim=dim,fieldName=:u)
    Φ=Field(materials,type=dim == 2 ? :ScalarField : :VectorField,dim=dim,fieldName=:φ)
    mpc_u=MPC(master="rightremote",slave="right",field=U)
    mpc_φ=MPC(master="rightremote",slave="right",field=Φ)
    R=rigidRotationMap(mpc_u,mpc_φ)
    Ku=∫(SymGrad(U)⋅D(dim == 2 ? :PlaneStress : :Solid,mat)⋅SymGrad(U);Ω="body")
    Mu=∫(U⋅U;Ω="body")
    Cu=0.2Mu
    block(A)=SystemMatrix([A A*R; R'*A R'*A*R])
    K,M,C=block(Ku),block(Mu),block(Cu)
    keys_u=dim == 2 ? (;ux=0.0,uy=0.0) : (;ux=0.0,uy=0.0,uz=0.0)
    support=[BoundaryCondition("left";field=U,keys_u...)]
    zu=LowLevelFEM._transient_field(U,zeros(ndofs(U),1),[0.0])
    zφ=LowLevelFEM._transient_field(Φ,zeros(ndofs(Φ),1),[0.0])
    initial=SystemVector([zu,zφ])
    mpc=[mpc_u,mpc_φ]
    prepared=LowLevelFEM.prepare_transient_system(K,M,C,initial,initial,initial,2,0.01;support,mpc)
    # Independent physical Galerkin system: u_phys = [I R] * T * q.
    physical_map=hcat(sparse(I,ndofs(U),ndofs(U)),R.A)*prepared.T
    G=physical_map[:,prepared.free]
    k,m,c=Matrix(G'*Ku.A*G),Matrix(G'*Mu.A*G),Matrix(G'*Cu.A*G)
    dt=0.1/sqrt(maximum(eigvals(Symmetric(k),Symmetric(m))))
    return (;U,Φ,R,Ku,Mu,Cu,K,M,C,support,mpc,initial,physical_map,G,k,m,c,dt)
end

function check_remote_kinematics(u,φ,fixture)
    (;U,Φ,mpc,R)=fixture
    for (field,P,condition) in ((u,U,mpc[1]),(φ,Φ,mpc[2]))
        rep=LowLevelFEM.mpcRepresentativeMap(P,[condition])
        tied=findall(i->rep[i]!=i,eachindex(rep))
        @test !isempty(tied)
        @test field.a[tied,:] ≈ field.a[rep[tied],:] atol=1e-11
    end
    master_tag=LowLevelFEM.getTagForPhysicalName("rightremote")
    nodes,coords=gmsh.model.mesh.getNodesForPhysicalGroup(-1,master_tag)
    master=Int(only(nodes))
    slave_tag=LowLevelFEM.getTagForPhysicalName("right")
    slaves,xyz=gmsh.model.mesh.getNodesForPhysicalGroup(-1,slave_tag)
    physical=u.a+R.A*φ.a
    for (i,node) in enumerate(slaves)
        lever=xyz[3i-2:3i]-coords
        urows=U.pdim*(master-1).+(1:U.pdim)
        rows=U.pdim*(Int(node)-1).+(1:U.pdim)
        if U.dim == 2
            expected=u.a[urows,:]+reshape([-lever[2],lever[1]],:,1)*φ.a[master:master,:]
        else
            φrows=Φ.pdim*(master-1).+(1:Φ.pdim)
            expected=hcat((u.a[urows,j]+cross(φ.a[φrows,j],lever) for j in axes(u.a,2))...)
        end
        @test physical[rows,:] ≈ expected atol=1e-11
    end
end

@testset "Transient remote translations and rotations (2D/3D)" begin
    for dim in (2,3)
        with_test_gmsh() do
            f=transient_remote_fixture(dim)
            (;U,Φ,R,Ku,Mu,Cu,K,M,C,support,mpc,initial,G,k,m,c,dt)=f
            @test isposdef(Symmetric(m))
            @test norm(R.A)>0
            n=21
            t=collect(0:n-1)*dt
            base=∫(U⋅(dim == 2 ? [0.0,1.0] : [0.0,1.0,0.2]);Γ="right")
            force=base.a*(1 .+ sin.(3t/dt/n))'
            fu=LowLevelFEM._transient_field(U,force,t)
            fφ=LowLevelFEM._transient_field(Φ,R.A'*force,t)
            F=SystemVector([fu,fφ])
            reduced_force=G'*force
            for run in (HHT,CDM)
                (u,φ),(v,ω)=run(K,M,C,F,initial,initial,n,dt;support,mpc)
                check_remote_kinematics(u,φ,f)
                physical=u.a+R.A*φ.a
                velocity=v.a+R.A*ω.a
                q=zeros(size(G,2)); qv=zero(q); a=m\reduced_force[:,1]
                previous=q-dt*qv+dt^2/2*a
                for j in 2:n
                    if run === HHT
                        qp=q+dt*qv+dt^2/4*a
                        vp=qv+dt/2*a
                        a=(m+dt/2*c+dt^2/4*k)\(reduced_force[:,j]-k*qp-c*vp)
                        q=qp+dt^2/4*a
                        qv=vp+dt/2*a
                    else
                        qnew=2q-previous+dt^2*(m\(reduced_force[:,j-1]-k*q-c*((q-previous)/dt)))
                        qv=(qnew-q)/dt
                        previous,q=q,qnew
                    end
                    @test physical[:,j] ≈ G*q rtol=1e-9 atol=1e-12
                    @test velocity[:,j] ≈ G*qv rtol=1e-9 atol=1e-12
                end
                @test !iszero(φ.a)
            end
            u,φ=FDM(K,M,F,initial,n,dt;support,mpc)
            check_remote_kinematics(u,φ,f)
            q=zeros(size(G,2))
            for j in 2:n
                q=(m+dt/2*k)\((m-dt/2*k)*q+dt/2*(reduced_force[:,j]+reduced_force[:,j-1]))
                @test u.a[:,j]+R.A*φ.a[:,j] ≈ G*q rtol=1e-9 atol=1e-12
            end
            # Exercise iterative solvers on the same remote block system.
            (ui,φi),_=HHT(Symmetric(K),Symmetric(M),Symmetric(C),F,initial,initial,n,dt;
                          support,mpc,solver=:cg,solveroptions=(reltol=1e-12,maxiter=300))
            (ud,φd),_=HHT(K,M,C,F,initial,initial,n,dt;support,mpc)
            @test ui.a+R.A*φi.a ≈ ud.a+R.A*φd.a rtol=1e-9 atol=1e-12
        end
    end
end

@testset "Remote rotation prescribed by a full U0 block history" begin
    for dim in (2,3)
        with_test_gmsh() do
            f=transient_remote_fixture(dim)
            (;U,Φ,R,Ku,Mu,Cu,K,M,C,mpc,initial,dt)=f
            n=15
            t=collect(0:n-1)*dt
            uk=dim == 2 ? (;ux=0.0,uy=0.0) : (;ux=0.0,uy=0.0,uz=0.0)
            φk=dim == 2 ? (;φ=0.0) : (;φx=0.0,φy=0.0,φz=0.0)
            support=vcat(f.support,[BoundaryCondition("rightremote";field=U,uk...),
                                    BoundaryCondition("rightremote";field=Φ,φk...)])
            remoteφ=constrainedDoFs(Φ,[support[end]])
            history=zeros(size(K.A,1),n)
            off=ndofs(U)
            history[off+remoteφ[end],:].=0.01t
            U0=SystemVector(history,nothing,K.problems,K.offsets)
            v0=zeros(size(K.A,1),1)
            v0[off.+remoteφ,1].=99.0
            V0=SystemVector(v0,nothing,K.problems,K.offsets)
            p=LowLevelFEM.prepare_transient_system(K,M,C,initial,U0,V0,n,dt;support,mpc)
            G=f.physical_map[:,p.free]
            @test !isempty(p.free)
            for run in (HHT,CDM)
                (u,φ),(v,ω)=run(K,M,C,initial,U0,V0,n,dt;support,mpc)
                @test φ.a[remoteφ[end],:] ≈ 0.01t atol=1e-13
                @test ω.a[remoteφ[end],:] ≈ fill(0.01,n) atol=1e-10
                check_remote_kinematics(u,φ,f)
                x=u.a+R.A*φ.a
                vel=v.a+R.A*ω.a
                if run === HHT
                    acc=zeros(size(x))
                    acc[:,1].=G*((G'*Mu.A*G)\(-G'*(Cu.A*vel[:,1]+Ku.A*x[:,1])))
                    for j in 2:n
                        acc[:,j].=(vel[:,j]-vel[:,j-1]-dt/2*acc[:,j-1])/(dt/2)
                        terms=(Mu.A*acc[:,j],Cu.A*vel[:,j],Ku.A*x[:,j])
                        scale=sum(norm(G'*term) for term in terms)
                        @test norm(G'*sum(terms)) < 1e-8*max(scale,1e-8)
                    end
                else
                    for j in 2:n-1
                        acc=(x[:,j+1]-2x[:,j]+x[:,j-1])/dt^2
                        terms=(Mu.A*acc,Cu.A*vel[:,j],Ku.A*x[:,j])
                        scale=sum(norm(G'*term) for term in terms)
                        @test norm(G'*sum(terms)) < 1e-8*max(scale,1e-8)
                    end
                end
            end
            u,φ=FDM(K,M,initial,U0,n,dt;support,mpc)
            @test φ.a[remoteφ[end],:] ≈ 0.01t atol=1e-13
            check_remote_kinematics(u,φ,f)
            x=u.a+R.A*φ.a
            for j in 2:n
                capacity=Mu.A*(x[:,j]-x[:,j-1])/dt
                stiffness=Ku.A*(x[:,j]+x[:,j-1])/2
                scale=norm(G'*capacity)+norm(G'*stiffness)
                @test norm(G'*(capacity+stiffness)) < 1e-8*max(scale,1e-8)
            end
        end
    end
end
