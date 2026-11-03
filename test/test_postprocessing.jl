@testset "Postprocessing: path time histories" begin
    with_test_gmsh() do
        structured_rect_mesh(nx=1, ny=1, order=2)
        times = [0.0, 0.2, 0.7]
        definitions = (
            (:ScalarField,
             (x,y,t) -> [1 + x + 2y + 3t],
             (x,y,t) -> 1 + x + 2y + 3t),
            (:VectorField,
             (x,y,t) -> [1 + x + t, 2 + y + 2t],
             (x,y,t) -> hypot(1 + x + t, 2 + y + 2t)),
            (:TensorField,
             (x,y,t) -> [1+x+t, 0.3+0.2t, 0, 0.3+0.2t, 2+y+t, 0, 0, 0, 0],
             (x,y,t) -> sqrt((1+x+t)^2 - (1+x+t)*(2+y+t) +
                              (2+y+t)^2 + 3*(0.3+0.2t)^2))
        )
        for (field_type, value, equivalent) in definitions
            P = Field([Material("body")]; type=field_type, dim=2, fieldName=:u)
            tags, coordinates, _ = gmsh.model.mesh.getNodes()
            values = zeros(ndofs(P), length(times))
            for (i, tag) in enumerate(tags), (j, time) in enumerate(times)
                rows = (Int(tag)-1)*P.pdim .+ (1:P.pdim)
                values[rows, j] .= value(coordinates[3i-2], coordinates[3i-1], time)
            end
            q = LowLevelFEM._transient_field(P, values, times)
            for field in (q, nodesToElements(q))
                source = isNodal(field) ? showDoFResults(field) : showElementResults(field)
                for step in (nothing, 0, 2)
                    tag, (x, y) = plotOnPath("bottom", source;
                        points=7, step, plot=true, offsetY=0.2)
                    types, counts, data = gmsh.view.getListData(tag)
                    selected = step === nothing ? eachindex(times) : (step+1:step+1)
                    samples = reshape(only(data), 3 + length(selected), 7)
                    expected = [equivalent(xi, 0.2, times[j])
                                for j in selected, xi in range(0, 1; length=7)]
                    @test types == ["SP"]
                    @test counts == [7]
                    @test samples[4:end, :] ≈ expected atol=1e-12
                    @test x ≈ collect(range(0, 1; length=7)) atol=1e-12
                    @test y ≈ vec(expected[1, :]) atol=1e-12
                    gmsh.view.remove(tag)
                end
                @test_throws ArgumentError plotOnPath("bottom", source; points=1)
                @test_throws ArgumentError plotOnPath("bottom", source; step=-1)
                @test_logs (:warn, r"using the last step") begin
                    tag = plotOnPath("bottom", source; points=2, step=99)
                    gmsh.view.remove(tag)
                end
                gmsh.view.remove(source)
            end
        end
    end
end
