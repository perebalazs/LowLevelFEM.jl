module LowLevelFEM

using LinearAlgebra, SparseArrays
using StaticArrays
using JLD2
using PrecompileTools: @setup_workload, @compile_workload
#using Base.Threads
import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh
export gmsh

include("debugtools.jl")
using .DebugTools
include("general.jl")
include("operators.jl")
include("linear.jl")
include("heat.jl")
include("nonlinear.jl")
include("nabla.jl")
include("poisson.jl")
include("multifield.jl")
include("constraints.jl")
include("solve.jl")
include("eigen.jl")
include("transient.jl")
include("transient_modify.jl")
include("spectral.jl")
include("contact.jl")
#include("fieldtools.jl")
#using .FieldTools
include("extra.jl")

export @showfields, @showstruct, @showdef, @showtype, @showmem, @showmethods, @disp, @showsize

# Representative public workflows, not an exhaustive API catalogue. Keep meshes
# tiny: the purpose is to capture JIT specializations, not to benchmark assembly.
# Setup and all native Gmsh state are confined to package precompilation.
const _PRECOMPILE_NEW = true
if _PRECOMPILE_NEW
    @setup_workload begin
        gmsh_was_initialized = gmsh.isInitialized() != 0
        gmsh_was_initialized || gmsh.initialize()
        terminal_pc = gmsh.option.getNumber("General.Terminal")
        verbosity_pc = gmsh.option.getNumber("General.Verbosity")

        try
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)

            # 3D elasticity and shared field/postprocessing kernels.
            structured_box_mesh(n=1, order=1)
            mat_pc = Material("body")
            @compile_workload begin
                U_pc = Field([mat_pc], type=:VectorField, dim=3, fieldName=:u)
                P_pc = Field([mat_pc], type=:ScalarField, dim=3, fieldName=:s)
                T_pc = Field([mat_pc], type=:TensorField, dim=3, fieldName=:T)

                K_pc = ∫(SymGrad(U_pc) ⋅ D(:Solid, mat_pc) ⋅ SymGrad(U_pc))
                ∫(U_pc ⋅ [1.0, 1.0, 1.0])
                f_pc = ∫(U_pc ⋅ [1.0, 0.0, 0.0]; Γ="right")
                support_pc = [BoundaryCondition("left"; field=U_pc, ux=0.0, uy=0.0, uz=0.0)]
                solveField(K_pc, f_pc; support=support_pc)
                solveField(Symmetric(K_pc), f_pc; support=support_pc)

                ∫(Curl(U_pc) ⋅ Curl(U_pc))
                ∫(TensorDiv(T_pc) ⋅ TensorDiv(T_pc))
                ∫(SurfaceGrad(P_pc) ⋅ SurfaceGrad(P_pc); Γ="right")

                s_pc = scalarField(P_pc, "body", 2.0)
                S_pc = ScalarField(P_pc, "body", 2.0)
                V_pc = VectorField(U_pc, "body", [1.0, 2.0, 3.0])
                F_pc = TensorField(U_pc, "body", [2.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0])
                ∇(s_pc)
                ∇(S_pc)
                ∇(V_pc)
                ∇(V_pc; nabla=:div)
                ∇(V_pc; nabla=:curl)
                ∇(F_pc; nabla=:div)
                ∂x(s_pc)
                ∂y(s_pc)
                ∂z(s_pc)
                nodesToElements(elementsToNodes(S_pc))
                elementsToElements(S_pc; onPhysicalGroup="left", fromPhysicalGroup="body")
                elementsToElements(V_pc; onPhysicalGroup="left")
                elementsToElements(F_pc; onPhysicalGroup="left")

                # Common coefficient and stress algebra, including the costly
                # small-matrix functions. Simple container predicates/constructors
                # need no separate dummy-data calls.
                S_pc + S_pc
                S_pc - S_pc
                S_pc * S_pc
                S_pc / S_pc
                sqrt(S_pc * S_pc + 1.0)
                log(S_pc)
                V_pc + V_pc
                V_pc - V_pc
                V_pc * S_pc
                V_pc / S_pc
                dot(V_pc, V_pc)
                V_pc × V_pc
                norm(V_pc)
                V_pc ∘ V_pc
                F_pc + F_pc
                F_pc - F_pc
                F_pc * F_pc
                F_pc ⋅ F_pc
                F_pc * V_pc
                V_pc * F_pc
                F_pc * S_pc
                F_pc / 2.0
                transpose(F_pc)
                unitTensor(F_pc)
                trace(F_pc)
                det(F_pc)
                inv(F_pc)
                eigen(F_pc)
                sqrt(F_pc)
                log(F_pc)
                ∫(P_pc, "body", S_pc)
                ∫(P_pc, "right", s_pc)
            end

            # Scalar/coefficient/compound DSL forms and mixed-field transformations.
            # P2 geometry also exercises the reduced P1 pressure reconstruction.
            structured_rect_mesh(n=2, order=2)
            @compile_workload begin
                mat_pc = Material("body")
                S_pc = Field([mat_pc], type=:ScalarField, dim=2, fieldName=:s)
                V_pc = Field([mat_pc], type=:VectorField, dim=2, fieldName=:v)
                P_pc = Field([mat_pc], type=:ScalarField, dim=2, fieldName=:p, reducedOrder=true)
                ∫(Grad(S_pc) ⋅ Grad(S_pc))
                ∫(S_pc ⋅ S_pc)
                ∫(S_pc ⋅ 1.0)
                c_pc = scalarField(S_pc, "body", 2.0)
                ∫(Grad(S_pc) ⋅ c_pc ⋅ Grad(S_pc))
                ∫(Grad(S_pc) ⋅ [c_pc 0.0; 0.0 c_pc] ⋅ Grad(S_pc))
                B_pc = Grad(S_pc) + [1.0; 0.0;;] ⋅ S_pc
                ∫(B_pc' ⋅ [2.0 0.0; 0.0 1.0] ⋅ B_pc)

                C_pc = D(:PlaneStress, mat_pc)
                K_pc = ∫(SymGrad(V_pc) ⋅ C_pc ⋅ SymGrad(V_pc); threads=1)
                ∫(SymGrad(V_pc) ⋅ C_pc ⋅ SymGrad(V_pc); assembly=:ijv)
                f_pc = ∫(V_pc ⋅ [1.0, 0.0])
                support_pc = [BoundaryCondition("left"; field=V_pc, vx=0.0, vy=0.0)]
                mpc_pc = MPC(master="rightbottom", slave="right", field=V_pc, vx=true, vy=true)
                solveField(K_pc, f_pc; support=support_pc, mpc=[mpc_pc])
                cs_pc = CoordinateSystem(VectorField(V_pc, "body", [0.0, 1.0, 0.0]))
                solveField(K_pc, Global(f_pc); support=support_pc, coordSys=[cs_pc])

                A_pc = ∫((SymGrad(V_pc) ⋅ SymGrad(V_pc)) * 2.0)
                B_pc = ∫(Div(V_pc) ⋅ P_pc)
                Z_pc = 0.0 * ∫(P_pc ⋅ P_pc)
                Kblock_pc = SystemMatrix([A_pc B_pc; B_pc' Z_pc])
                Fblock_pc = SystemVector([f_pc, ∫(P_pc ⋅ 0.0)])
                support_block_pc = [
                    BoundaryCondition("top"; field=V_pc, vx=0.0, vy=0.0),
                    BoundaryCondition("bottom"; field=V_pc, vx=0.0, vy=0.0),
                    BoundaryCondition("rightbottom"; field=P_pc, p=0.0)
                ]
                solveField(Kblock_pc, Fblock_pc; support=support_block_pc)
            end

            # Generalized eigenproblem, using the built-in direct backend. Optional
            # ARPACK/iterative/AD extensions are not dependencies of this workload.
            line_mesh(n=2, order=1)
            @compile_workload begin
                P_pc = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
                K_pc = ∫(Grad(P_pc) ⋅ Grad(P_pc))
                M_pc = ∫(P_pc ⋅ P_pc)
                ∫(AxialGrad(P_pc) ⋅ AxialGrad(P_pc))
                support_pc = [BoundaryCondition("left"; field=P_pc, p=0.0),
                            BoundaryCondition("right"; field=P_pc, p=0.0)]
                modes_pc = solveEigen(K_pc, M_pc; n=1, solver=:direct, support=support_pc)
                solveEigen(Symmetric(K_pc), Symmetric(M_pc); n=1, solver=:direct,
                        support=support_pc, normalize=:mass)
                eigenValue(modes_pc)
                eigenVector(modes_pc)
            end

            # Two one-element patches exercise closest-point search, frozen-geometry
            # penalty assembly, L2 gap projection and multiplier coupling.
            gmsh.clear()
            gmsh.model.add("precompile_contact")
            lower_pc = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
            upper_pc = gmsh.model.occ.addRectangle(0.0, 1.1, 0.0, 1.0, 1.0)
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup(2, [lower_pc, upper_pc], -1, "body")
            for (name_pc, patch_pc, y_pc) in (("master", lower_pc, 1.0), ("slave", upper_pc, 1.1))
                sides_pc = gmsh.model.getBoundary([(2, patch_pc)], false, false)
                tags_pc = [tag for (d, tag) in sides_pc
                        if abs(gmsh.model.occ.getCenterOfMass(d, tag)[2] - y_pc) < 1e-10]
                gmsh.model.addPhysicalGroup(1, tags_pc, -1, name_pc)
            end
            for (_, tag_pc) in gmsh.model.getEntities(1)
                gmsh.model.mesh.setTransfiniteCurve(tag_pc, 2)
            end
            for (_, tag_pc) in gmsh.model.getEntities(2)
                gmsh.model.mesh.setTransfiniteSurface(tag_pc)
                gmsh.model.mesh.setRecombine(2, tag_pc)
            end
            gmsh.option.setNumber("Mesh.ElementOrder", 1)
            gmsh.model.mesh.generate(2)
            @compile_workload begin
                U_pc = Field([Material("body")], type=:VectorField, dim=2, fieldName=:u)
                u_pc = vectorField(U_pc, "body", [0.0, 0.0, 0.0])
                contact_pc = contact(U_pc; master="master", slave="slave", displacement=u_pc)
                G_pc = ContactGap(contact_pc; active=:all)
                ∫(G_pc ⋅ 2.0 ⋅ G_pc)
                Gfull_pc = ContactGap(contact_pc; components=:all, active=:all)
                ∫(Gfull_pc ⋅ ContactStiffness(contact_pc; cn=2.0, ct=0.0) ⋅ Gfull_pc)
                ContactGap(contact_pc, u_pc; gauss=:full)
                updateContact!(contact_pc, u_pc)
                Λ_pc = Field([Material("body")], type=:ScalarField, dim=2, fieldName=:λ)
                ∫(Λ_pc ⋅ G_pc)
            end
        finally
            if gmsh.isInitialized() != 0
                gmsh.clear()
                if gmsh_was_initialized
                    gmsh.option.setNumber("General.Terminal", terminal_pc)
                    gmsh.option.setNumber("General.Verbosity", verbosity_pc)
                else
                    gmsh.finalize()
                end
            end
        end
    end
end

# Legacy compatibility is intentionally isolated. Set this to false to omit
# only this block on the next package precompilation; it adds no runtime work.
const _PRECOMPILE_LEGACY = true
if _PRECOMPILE_LEGACY
    @setup_workload begin
        gmsh_was_initialized = gmsh.isInitialized() != 0
        gmsh_was_initialized || gmsh.initialize()
        terminal_pc = gmsh.option.getNumber("General.Terminal")
        verbosity_pc = gmsh.option.getNumber("General.Verbosity")
        try
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)
            structured_rect_mesh(n=1, order=1)
            @compile_workload begin
                Plegacy_pc = Problem([Material("body")], type=:PlaneStress)
                support_legacy_pc = [displacementConstraint("left"; ux=0.0),
                                     displacementConstraint("bottom"; uy=0.0)]
                u_legacy_pc = solveDisplacement(Plegacy_pc;
                    load=[load("right"; fx=1.0)], support=support_legacy_pc)
                solveStress(u_legacy_pc)
            end
        finally
            if gmsh.isInitialized() != 0
                gmsh.clear()
                if gmsh_was_initialized
                    gmsh.option.setNumber("General.Terminal", terminal_pc)
                    gmsh.option.setNumber("General.Verbosity", verbosity_pc)
                else
                    gmsh.finalize()
                end
            end
        end
    end
end

end #module
