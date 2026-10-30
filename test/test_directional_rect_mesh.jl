@testset "Structured rectangle: independent direction counts" begin
    with_test_gmsh() do
        structured_rect_mesh(lx=1.0,ly=0.04,nx=25,ny=3,order=2)
        P=Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        @test P.non == 51*7
        @test ndofs(P) == 714
        @test length(nodesToElements(initialDisplacement(P,"body";ux=0.0,uy=0.0)).numElem) == 75
        structured_rect_mesh(lx=1.0,ly=0.04,n=3,order=2)
        P=Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        @test P.non == 151*7
        @test ndofs(P) == 2114
        structured_rect_mesh(lx=1.0,ly=0.04,n=3,nx=25,order=2)
        P=Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        @test P.non == 51*7
        structured_rect_mesh(lx=1.0,ly=0.04,n=3,ny=1,order=2)
        P=Field([Material("body")],type=:VectorField,dim=2,fieldName=:u)
        @test P.non == 151*3
        structured_rect_mesh(nx=2,ny=1,leftRemote=true,rightRemote=true)
        @test length(gmsh.model.getEntitiesForPhysicalName("leftremote")) == 1
        @test length(gmsh.model.getEntitiesForPhysicalName("rightremote")) == 1
        for invalid in (0,-1,2.5,true)
            @test_throws ArgumentError structured_rect_mesh(nx=invalid)
            @test_throws ArgumentError structured_rect_mesh(ny=invalid)
        end
    end
end
