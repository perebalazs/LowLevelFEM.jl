@testset "Differential and tensor operators" begin
    with_test_gmsh() do
        structured_box_mesh(n=1, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=3, fieldName=:p)

        @testset "Constant and affine spatial derivatives" begin
            c = ScalarField(P, "body", 3.0)
            @test isElementwise(grad(c))
            test_constant_field(grad(c), zeros(3))

            affine(x, y, z) = 1.0 + 2x - 3y + 4z
            for f in (scalarField(P, "body", affine), ScalarField(P, "body", affine))
                test_constant_field(∇(f), [2.0, -3.0, 4.0])
                test_constant_field(∂x(f), 2.0)
                test_constant_field(∂y(f), -3.0)
                test_constant_field(∂z(f), 4.0)
            end

            # A nonsymmetric gradient catches swapped derivative/component axes.
            v = VectorField(P, "body", [(x,y,z)->2x+3y,
                                        (x,y,z)->-x+4z,
                                        (x,y,z)->5y-z])
            test_constant_field(∇(v), [2.0 3.0 0.0; -1.0 0.0 4.0; 0.0 5.0 -1.0])
            test_constant_field(∇(v; nabla=:div), 1.0)
            test_constant_field(∇(v; nabla=:curl), [1.0, 0.0, -4.0])

            T = TensorField(P, "body", [((x,y,z)->2x) ((x,y,z)->3y) 0.0;
                                         0.0 ((x,y,z)->-y) ((x,y,z)->4z);
                                         ((x,y,z)->5x) 0.0 ((x,y,z)->6z)])
            test_constant_field(∇(T; nabla=:div), [5.0, 3.0, 11.0])
        end

        @testset "Exact tensor and vector algebra" begin
            v = VectorField(P, "body", [1.0, 2.0, 3.0])
            w = VectorField(P, "body", [3.0, -1.0, 2.0])
            T0 = [2.0 1.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0]
            T = TensorField(P, "body", T0)
            test_constant_field(dot(v, w), 7.0)
            test_constant_field(v × w, [7.0, 7.0, -7.0])
            test_constant_field(norm(v), sqrt(14.0))
            test_constant_field(T * v, T0 * [1.0, 2.0, 3.0])
            test_constant_field(transpose(T), T0')
            test_constant_field(trace(T), 9.0)
            test_constant_field(det(T), 24.0)
            test_constant_field(inv(T), inv(T0))
            test_constant_field(T * inv(T), Matrix{Float64}(I, 3, 3))

            H = TensorField(P, "body", [4.0 0.0 0.0; 0.0 9.0 0.0; 0.0 0.0 16.0])
            test_constant_field(sqrt(H), [2.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0])
            values, vectors = eigen(H)
            # Tensor-field eigenvalues are returned in descending order.
            test_constant_field(values, [16.0, 9.0, 4.0])
            test_constant_field(vectors * transpose(vectors), Matrix{Float64}(I, 3, 3))
        end
    end
end

@testset "DSL operator conventions and energies" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1, order=1)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=2, fieldName=:u)
        # grad(u) = [1 2; 3 4], engineering strain = [1, 4, 5].
        u = vectorField(U, "body", [(x,y,z)->x+2y, (x,y,z)->3x+4y, 0.0])
        for (operator, energy) in ((Grad(U), 30.0), (SymGrad(U), 42.0),
                                   (Div(U), 25.0), (Curl(U), 1.0))
            K = ∫(operator ⋅ operator)
            @test dot(u.a[:, 1], K.A * u.a[:, 1]) ≈ energy atol=1e-10
        end
        # A diagonal coefficient distinguishes the component-major ordering
        # [ux,x, ux,y, uy,x, uy,y] from other flattening conventions.
        K = ∫(Grad(U) ⋅ Diagonal([1.0, 2.0, 3.0, 4.0]) ⋅ Grad(U))
        @test dot(u.a[:, 1], K.A * u.a[:, 1]) ≈ 100.0 atol=1e-10

        S = Field([mat], type=:ScalarField, dim=2, fieldName=:s)
        test_constant_field(grad(ScalarField(S, "body", 3.0)), zeros(3))
        test_constant_field(∇(scalarField(S, "body", (x,y,z)->2x-3y)), [2.0, -3.0, 0.0])

        T = Field([mat], type=:TensorField, dim=2, fieldName=:T)
        σ = tensorField(T, "body", [((x,y,z)->x) 0.0 0.0;
                                     0.0 ((x,y,z)->2y) 0.0; 0.0 0.0 0.0])
        # The full nine-component storage is retained for a 2D tensor field.
        Kdiv = ∫(TensorDiv(T) ⋅ TensorDiv(T))
        @test dot(σ.a[:, 1], Kdiv.A * σ.a[:, 1]) ≈ 5.0 atol=1e-10

        # Nonsymmetry distinguishes div(σ) from div(σ') in 2D as well.
        σ = tensorField(T, "body", [((x,y,z)->2x) ((x,y,z)->3y) 0.0;
                                     0.0 ((x,y,z)->-y) 0.0; 0.0 0.0 0.0])
        test_constant_field(∇(σ; nabla=:div), [5.0, -1.0, 0.0])
        for assembly in (:csc, :ijv), threads in (1, :auto)
            Kdiv = ∫(TensorDiv(T) ⋅ TensorDiv(T); assembly=assembly, threads=threads)
            @test dot(σ.a[:, 1], Kdiv.A * σ.a[:, 1]) ≈ 26.0 atol=1e-10
        end
    end

    with_test_gmsh() do
        structured_box_mesh(n=1, order=1)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=3, fieldName=:u)
        u = vectorField(U, "body", [(x,y,z)->y+2z, (x,y,z)->3z, 0.0])
        # Documented 3D engineering order: [xx, yy, zz, xy, yz, zx].
        strain = [0.0, 0.0, 0.0, 1.0, 3.0, 2.0]
        C = Diagonal([1.0, 2.0, 3.0, 5.0, 7.0, 11.0])
        K = ∫(SymGrad(U) ⋅ C ⋅ SymGrad(U))
        # Unequal shear coefficients expose an incorrect Voigt permutation.
        @test dot(u.a[:, 1], K.A * u.a[:, 1]) ≈ dot(strain, C * strain) atol=1e-10

        T = Field([mat], type=:TensorField, dim=3, fieldName=:T)
        σ = tensorField(T, "body", [((x,y,z)->2x) ((x,y,z)->3y) 0.0;
                                     0.0 ((x,y,z)->-y) ((x,y,z)->4z);
                                     ((x,y,z)->5x) 0.0 ((x,y,z)->6z)])
        # div(σ)=[5,3,11], with the documented column-major tensor storage.
        for assembly in (:csc, :ijv), threads in (1, :auto)
            Kdiv = ∫(TensorDiv(T) ⋅ TensorDiv(T); assembly=assembly, threads=threads)
            @test dot(σ.a[:, 1], Kdiv.A * σ.a[:, 1]) ≈ 155.0 atol=1e-10
        end
    end
end
