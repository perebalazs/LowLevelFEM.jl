var documenterSearchIndex = {"docs":
[{"location":"Examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"Examples/#2D-Cantilever","page":"Examples","title":"2D Cantilever","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: $\\sigma_x$ on deformed shape)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: $\\sigma_x$ and $\\tau_{yx}$ on path)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever2D.jl","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"import LowLevelFEM as FEM\nusing LowLevelFEM\n\ngmsh.initialize()\n\ngmsh.open(\"cantilever2D.geo\")\nmat = FEM.material(\"body\", E=2.e5, ν=0.3)\nproblem = FEM.Problem([mat], type=\"PlaneStress\")\n\nsupp = FEM.displacementConstraint(\"supp\", ux=0, uy=0)\nload = FEM.load(\"load\", fy=-1)\n\nK = FEM.stiffnessMatrix(problem)\nf = FEM.loadVector(problem, [load])\n\nFEM.applyBoundaryConditions!(problem, K, f, [supp])\n\nq = FEM.solveDisplacement(K, f)\nS = FEM.solveStress(problem, q)\n\nu = FEM.showDoFResults(problem, q, \"uvec\", name=\"uvec\", visible=false)\nux = FEM.showDoFResults(problem, q, \"ux\", name=\"ux\", visible=false)\nuy = FEM.showDoFResults(problem, q, \"uy\", name=\"uy\", visible=false)\n\ns = FEM.showStressResults(problem, S, \"s\", name=\"σ\", visible=true, smooth=true)\nsx = FEM.showStressResults(problem, S, \"sx\", name=\"σx\", visible=false, smooth=true)\nsy = FEM.showStressResults(problem, S, \"sy\", name=\"σy\", visible=false, smooth=true)\nsxy = FEM.showStressResults(problem, S, \"sxy\", name=\"τxy\", visible=false, smooth=true)\n\nFEM.plotOnPath(problem, \"path\", sx, 100, name=\"σx\", visible=false);\nFEM.plotOnPath(problem, \"path\", sxy, 100, name=\"τxy\", visible=false);\nFEM.plotOnPath(problem, \"path\", ux, 100, name=\"ux\", visible=false);\n\ngmsh.fltk.run()\ngmsh.finalize()","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever2D.geo","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"SetFactory(\"OpenCASCADE\");\n\nRectangle(1) = {0, 0, 0, 100, 10, 0};\n\nPhysical Curve(\"supp\", 5) = {4};\nPhysical Curve(\"load\", 6) = {2};\nPhysical Surface(\"body\", 7) = {1};\n\nRecombine Surface {1};\n\nTransfinite Line {2,4} = 4;\nTransfinite Line {1,3} = 31;\nTransfinite Surface {1};\n\nMesh.ElementOrder = 3;\n\nSetName \"cantilever2D\";\nMesh 2;\n\nPoint(5) = {10, 0, 0, 1.0};\nPoint(6) = {10, 10, 0, 1.0};\nLine(5) = {5, 6};\n\nPhysical Curve(\"path\", 8) = {5};","category":"page"},{"location":"Examples/#3D-Cantilever","page":"Examples","title":"3D Cantilever","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: $\\sigma_x$ on deformed shape)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever3D.jl","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"import LowLevelFEM as FEM\nusing LowLevelFEM\n\ngmsh.initialize()\n\ngmsh.open(\"cantilever3D.geo\")\nmat = FEM.material(\"body\", E=2.e5, ν=0.3)\nproblem = FEM.Problem([mat])\n\nsupp = FEM.displacementConstraint(\"supp\", ux=0, uy=0, uz=0)\nload = FEM.load(\"load\", fy=-1)\n\nK = FEM.stiffnessMatrix(problem)\nf = FEM.loadVector(problem, [load])\n\nFEM.applyBoundaryConditions!(problem, K, f, [supp])\n\nq = FEM.solveDisplacement(K, f)\nS = FEM.solveStress(problem, q)\n\nu = FEM.showDoFResults(problem, q, \"uvec\", name=\"uvec\", visible=false)\nux = FEM.showDoFResults(problem, q, \"ux\", name=\"ux\", visible=false)\nuy = FEM.showDoFResults(problem, q, \"uy\", name=\"uy\", visible=false)\nuz = FEM.showDoFResults(problem, q, \"uz\", name=\"uz\", visible=false)\n\ns = FEM.showStressResults(problem, S, \"s\", name=\"σ\", visible=true, smooth=true)\nsx = FEM.showStressResults(problem, S, \"sx\", name=\"σx\", visible=false, smooth=true)\nsy = FEM.showStressResults(problem, S, \"sy\", name=\"σy\", visible=false, smooth=true)\nsz = FEM.showStressResults(problem, S, \"sz\", name=\"σz\", visible=false, smooth=true)\nsxy = FEM.showStressResults(problem, S, \"sxy\", name=\"τxy\", visible=false, smooth=true)\nsyz = FEM.showStressResults(problem, S, \"syz\", name=\"τyz\", visible=false, smooth=true)\nszx = FEM.showStressResults(problem, S, \"szx\", name=\"τzx\", visible=false, smooth=true)\n\nFEM.plotOnPath(problem, \"path\", sx, 100, name=\"σx\", visible=false);\nFEM.plotOnPath(problem, \"path\", sxy, 100, name=\"τxy\", visible=false);\nFEM.plotOnPath(problem, \"path\", ux, 100, name=\"ux\", visible=false);\n\ngmsh.fltk.run()\ngmsh.finalize()","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever3D.geo","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"SetFactory(\"OpenCASCADE\");\n\nBox(1) = {0, 0, 0, 100, 10, 10};\n\nPhysical Surface(\"supp\", 13) = {1};\nPhysical Surface(\"load\", 14) = {2};\nPhysical Volume(\"body\", 15) = {1};\n\nRecombine Surface {1:6};\n\nTransfinite Line {1:8} = 4;\nTransfinite Line {9:12} = 31;\nTransfinite Surface {1:6};\nTransfinite Volume {1};\n\nMesh.ElementOrder = 3;\n\nSetName \"cantilever3D\";\nMesh 3;\n\nPoint(9) = {10, 0, 5, 1.0};\nPoint(10) = {10, 10, 5, 1.0};\nLine(13) = {9, 10};\n\nPhysical Curve(\"path\", 16) = {13};","category":"page"},{"location":"Examples/#L-shaped-plate","page":"Examples","title":"L-shaped plate","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: Mesh with a path for graphs)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: Fillet)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: Equivalent stress)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: Equivalent stress on path)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"LshapedPlate.jl","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"import LowLevelFEM as FEM\nusing LowLevelFEM\n\ngmsh.initialize()\n\n#gmsh.open(\"LshapedPlate.geo\")\ngmsh.open(\"LshapedPlate2.geo\")\n\nmat = FEM.material(\"body\", E=2.e5, ν=0.3)\nproblem = FEM.Problem([mat], type=\"PlaneStress\", thickness=1)\n\nbc1 = FEM.displacementConstraint(\"fix\", ux=0, uy=0)\nld1 = FEM.load(\"load\", fy=-1)\n\nK = FEM.stiffnessMatrix(problem)\nf = FEM.loadVector(problem, [ld1])\nFEM.applyBoundaryConditions!(problem, K, f, [bc1])\n\nq = FEM.solveDisplacement(K, f)\nS = FEM.solveStress(problem, q)\n\nu = FEM.showDoFResults(problem, q, \"uvec\", name=\"uvec\", visible=false)\nux = FEM.showDoFResults(problem, q, \"ux\", name=\"ux\", visible=false)\nuy = FEM.showDoFResults(problem, q, \"uy\", name=\"uy\", visible=false)\nuz = FEM.showDoFResults(problem, q, \"uz\", name=\"uz\", visible=false)\ns = FEM.showStressResults(problem, S, \"s\", name=\"σ red\", visible=false, smooth=false)\nss = FEM.showStressResults(problem, S, \"s\", name=\"σ red smooth\", visible=true, smooth=true)\nsx = FEM.showStressResults(problem, S, \"sx\", name=\"σx\", visible=false, smooth=true)\nsy = FEM.showStressResults(problem, S, \"sy\", name=\"σy\", visible=false, smooth=true)\nsz = FEM.showStressResults(problem, S, \"sz\", name=\"σz\", visible=false, smooth=true)\nsxy = FEM.showStressResults(problem, S, \"sxy\", name=\"τxy\", visible=false, smooth=true)\nsyz = FEM.showStressResults(problem, S, \"syz\", name=\"τyz\", visible=false, smooth=true)\nszx = FEM.showStressResults(problem, S, \"szx\", name=\"τzx\", visible=false, smooth=true)\n\nFEM.plotOnPath(problem, \"path\", s, 100, name=\"σred\", visible=false);\n\ngmsh.fltk.run()\ngmsh.finalize()","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"LshapedPlate.geo","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Point(1) = {0, 0, 0, 15.0};\nPoint(2) = {100, 0, 0, 15.0};\nPoint(3) = {100, 50, 0, 15.0};\nPoint(4) = {50, 50, 0, 0.5};\nPoint(5) = {50, 100, 0, 15.0};\nPoint(6) = {0, 100, 0, 15.0};\nLine(1) = {1, 2};\nLine(2) = {2, 3};\nLine(3) = {3, 4};\nLine(4) = {4, 5};\nLine(5) = {5, 6};\nLine(6) = {6, 1};\nCurve Loop(1) = {6, 1, 2, 3, 4, 5};\nPlane Surface(1) = {1};\n\nPhysical Curve(\"fix\", 7) = {5};\nPhysical Curve(\"load\", 8) = {2};\nPhysical Surface(\"body\", 11) = {1};\n\nSetName \"Lshape\";\n\nMesh.ElementOrder = 4;\nMesh.HighOrderOptimize = 1;\nMesh 2;\n\nPoint(7) = {0, 0, 0, 1.0};\nPoint(8) = {50, 50, 0, 1.0};\nLine(7) = {7, 8};\n\nPhysical Curve(\"path\", 9) = {7};","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"LshapedPlate2.geo","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"R=1;\n\nPoint(1) = {0, 0, 0, 15.0};\nPoint(2) = {100, 0, 0, 15.0};\nPoint(3) = {100, 50, 0, 15.0};\nPoint(4) = {50+R, 50, 0, R/1.6};\nPoint(5) = {50, 50+R, 0, R/1.6};\nPoint(6) = {50, 100, 0, 15.0};\nPoint(7) = {0, 100, 0, 15.0};\nPoint(8) = {50+R, 50+R, 0, 0.0};\nLine(1) = {1, 2};\nLine(2) = {2, 3};\nLine(3) = {3, 4};\nCircle(4) = {4, 8, 5};\nLine(5) = {5, 6};\nLine(6) = {6, 7};\nLine(7) = {7, 1};\nCurve Loop(1) = {1, 2, 3, 4, 5, 6, 7};\nPlane Surface(1) = {1};\n\nPhysical Curve(\"fix\", 8) = {6};\nPhysical Curve(\"load\", 9) = {2};\nPhysical Surface(\"body\", 11) = {1};\n\nSetName \"Lshape\";\nMesh.ElementOrder = 4;\nMesh.HighOrderOptimize = 1;\nMesh 2;\n\nPoint(9) = {0, 0, 0, 1.0};\nPoint(10) = {50+0.415*R, 50+0.415*R, 0, 1.0};\nLine(8) = {9, 10};\n\nPhysical Curve(\"path\", 10) = {8};","category":"page"},{"location":"Examples/#Wave-propagation-in-a-plate","page":"Examples","title":"Wave propagation in a plate","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: velocity field)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"wavePropagation.jl","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"import LowLevelFEM as FEM\nusing LowLevelFEM\n\ngmsh.initialize()\n\nE = 2e5\nν = 0.3\nρ = 7.85e-9\nthick = 1\nheight = 10\nbase = 100\nelemSize = 2 #22\n\napproxOrder = 2\ninternalNodes = true\nquadElements = true\n\ngmsh.model.add(\"rectangle\")\n\np1 = gmsh.model.occ.addPoint(0, 0, 0)\np2 = gmsh.model.occ.addPoint(base, 0, 0)\np3 = gmsh.model.occ.addPoint(base, height, 0)\np4 = gmsh.model.occ.addPoint(0, height, 0)\n\nl1 = gmsh.model.occ.addLine(p1, p2)\nl2 = gmsh.model.occ.addLine(p2, p3)\nl3 = gmsh.model.occ.addLine(p3, p4)\nl4 = gmsh.model.occ.addLine(p4, p1)\n\ncl1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])\n\nl5 = gmsh.model.occ.addCircle(base / 2, height / 2, 0, min(base, height) / 4)\ncl2 = gmsh.model.occ.addCurveLoop([l5])\n\nsf1 = gmsh.model.occ.addPlaneSurface([cl1, cl2])\n\ngmsh.model.occ.synchronize()\n\nphg = gmsh.model.addPhysicalGroup(1, [l2])\ngmsh.model.setPhysicalName(1, phg, \"supp\")\nphg = gmsh.model.addPhysicalGroup(1, [l4])\ngmsh.model.setPhysicalName(1, phg, \"load\")\nphg = gmsh.model.addPhysicalGroup(2, [sf1])\ngmsh.model.setPhysicalName(2, phg, \"body\")\n\nFEM.generateMesh(sf1, elemSize, approxOrder=approxOrder, algorithm=6, quadrangle=quadElements, internalNodes=internalNodes)\n\nmat = FEM.material(\"body\", E=E, ν=ν)\nproblem = FEM.Problem([mat], type=\"PlaneStress\", thickness=thick)\n\nsupp = FEM.displacementConstraint(\"supp\", ux=0, uy=0)\nload = FEM.load(\"load\", fx=1, fy=0)\n\ngmsh.option.setNumber(\"Mesh.Lines\", 0)\n\nK = FEM.stiffnessMatrix(problem)\nf = FEM.loadVector(problem, [load])\nM = FEM.massMatrix(problem)\nC = 4e-3 * K\n\nFEM.applyBoundaryConditions!(problem, K, M, C, f, [supp]);\n\nTₘᵢₙ = FEM.smallestPeriodTime(K, M)\nq = FEM.solveDisplacement(K, f)\n\ndof, dof = size(K)\nu0 = zeros(dof)\nv0 = zeros(dof)\nFEM.initialDisplacement!(problem, \"supp\", u0, ux=0)\nFEM.initialVelocity!(problem, \"body\", v0, vx=1000)\nFEM.initialVelocity!(problem, \"supp\", v0, vx=0)\nf = zeros(dof)\n\nE = problem.material[1][2]\nρ = problem.material[1][4]\nc = √(E / ρ)\nξₘₐₓ = 1e-1\nβ = ξₘₐₓ * Tₘᵢₙ / π\nC = β * K\nu, v, t = FEM.CDM(K, M, C, f, u0, v0, base / c * 2, Tₘᵢₙ / π * (√(1 + ξₘₐₓ^2) - ξₘₐₓ) * 1.0)\n\nS = FEM.solveStress(problem, q)\n\nuvec = FEM.showDoFResults(problem, q, \"uvec\", name=\"uvec\", visible=false)\nux = FEM.showDoFResults(problem, q, \"ux\", name=\"ux\", visible=false)\nuy = FEM.showDoFResults(problem, q, \"uy\", name=\"uy\", visible=false)\nuz = FEM.showDoFResults(problem, q, \"uz\", name=\"uz\", visible=false)\ns = FEM.showStressResults(problem, S, \"s\", name=\"σ\", visible=false, smooth=true)\nsx = FEM.showStressResults(problem, S, \"sx\", name=\"σx\", visible=false, smooth=true)\nsy = FEM.showStressResults(problem, S, \"sy\", name=\"σy\", visible=false, smooth=true)\nsz = FEM.showStressResults(problem, S, \"sz\", name=\"σz\", visible=false, smooth=true)\nsxy = FEM.showStressResults(problem, S, \"sxy\", name=\"τxy\", visible=false, smooth=true)\nsyz = FEM.showStressResults(problem, S, \"syz\", name=\"τyz\", visible=false, smooth=true)\nszx = FEM.showStressResults(problem, S, \"szx\", name=\"τzx\", visible=false, smooth=true)\nvvec = FEM.showDoFResults(problem, v, t=t, \"uvec\", name=\"v(t)\", visible=true)\ngmsh.view.option.setNumber(vvec, \"NormalRaise\", 0.03)\n\nsts = ceil(Int64, (base / c * 2) / 6 / (Tₘᵢₙ / π * (√(1 + ξₘₐₓ^2) - ξₘₐₓ)))\ndisplay(sts)\nSp = FEM.solveStress(problem, u[:, sts])\nsp = FEM.showStressResults(problem, Sp, \"s\", name=\"σ at t\", visible=false, smooth=false);\n\nSanim = FEM.solveStress(problem, u[:, 1:sts])\nsanim = FEM.showStressResults(problem, Sanim, \"s\", t=t[1:sts], name=\"σ anim\", visible=false, smooth=false);\n\ngmsh.fltk.run()\ngmsh.finalize()","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"For more examples see examples on GitHub","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"Functions/#LowLevelFEM.jl","page":"Functions","title":"LowLevelFEM.jl","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Documentation for LowLevelFEM.jl","category":"page"},{"location":"Functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Modules = [LowLevelFEM]","category":"page"},{"location":"Functions/#LowLevelFEM.Problem","page":"Functions","title":"LowLevelFEM.Problem","text":"Problem(names; thickness=..., type=..., bandwidth=...)\n\nA structure containing the most important data of the problem. \n\nname of the model (in gmsh)\ntype of the problem: 3D \"Solid\", \"PlaneStrain\" or \"PlaneStress\"\nbandwidth optimization using built-in gmsh function. Possibilities: \"RCMK\" (default), Hilbert\" and \"Metis\"\ndimension of the problem, determined from type\nmaterial constants: Physical group, Young's modulus, Poisson's ratio, mass density (in vector of tuples names)\nthickness of the plate\nnumber of nodes (non)\n\nTypes:\n\nnames: Vector{Touple{String, Float64, Float64, Float64}}\ntype: String\nbandwidth: String\ndim: Integer\nthickness: Float64\nnon: Integer\n\n\n\n\n\n","category":"type"},{"location":"Functions/#LowLevelFEM.StressField","page":"Functions","title":"LowLevelFEM.StressField","text":"StressField(sigma, numElem, nsteps)\n\nA structure containing the data of a stress field. \n\nsigma: vector of ElementNodeData type stress data (see gmsh.jl)\nnumElem: vector of tags of elements\nnsteps: number of stress fields stored in sigma (for animations).\n\nTypes:\n\nsigma: Vector{Matrix{Float64}}\nnumElem: Vector{Integer}\nnsteps: Integer\n\n\n\n\n\n","category":"type"},{"location":"Functions/#LowLevelFEM.CDM-NTuple{8, Any}","page":"Functions","title":"LowLevelFEM.CDM","text":"FEM.CDM(K, M, C, f, u0, v0, T, Δt)\n\nSolves a transient dynamic problem using central difference method (explicit). K is the stiffness Matrix, M is the mass matrix, C is the damping matrix, f is the load vector, u0 is the initial displacement, v0 is the initial velocity, T is the upper bound ot the time intervall (lower bound is zero) and Δt is the time step size. Returns the displacement vectors and velocity vectors in each time step arranged in the columns of the two matrices u and v and a vector t of the time instants used.\n\nReturn: u, v, t\n\nTypes:\n\nK: SparseMatrix\nM: SparseMatrix\nC: SparseMatrix\nf: Vector{Float64}\nu0: Vector{Float64}\nv0: Vector{Float64}\nT: Float64\nΔt: Float64 \nu: Matrix{Float64}\nv: Matrix{Float64}\nt: Vector{Float64}\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.applyBoundaryConditions!-NTuple{4, Any}","page":"Functions","title":"LowLevelFEM.applyBoundaryConditions!","text":"FEM.applyBoundaryConditions!(problem, stiffMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nstiffMat: SparseMatrix \nloadVec: Vector \nsupports: Vector{Tuple{String, Float64, Float64, Float64}}\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.applyBoundaryConditions!-NTuple{6, Any}","page":"Functions","title":"LowLevelFEM.applyBoundaryConditions!","text":"FEM.applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat, mass matrix massMat, damping matrix dampMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nstiffMat: SparseMatrix \nmassMat: SparseMatrix \ndampMat: SparseMatrix \nloadVec: Vector{Float64}\nsupports: Vector{Tuple{String, Float64, Float64, Float64}}\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.applyBoundaryConditions-NTuple{4, Any}","page":"Functions","title":"LowLevelFEM.applyBoundaryConditions","text":"FEM.applyBoundaryConditions(problem, stiffMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz. Creates a new stiffness matrix and load vector.\n\nReturn: stiffMat, loadVec\n\nTypes:\n\nproblem: Problem\nstiffMat: SparseMatrix \nloadVec: Vector \nsupports: Vector{Tuple{String, Float64, Float64, Float64}}\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.displacementConstraint-Tuple{Any}","page":"Functions","title":"LowLevelFEM.displacementConstraint","text":"FEM.displacementConstraint(name; ux=..., uy=..., uz=...)\n\nGives the displacement constraints on name physical group. At least one ux,  uy or uz value have to be given (depending on the dimension of the problem).\n\nReturn: none\n\nTypes:\n\nname: string\nux: double\nuy: double\nuz: double\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.generateMesh-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.generateMesh","text":"FEM.generateMesh(problem, surf, elemSize; approxOrder=1, algorithm=6, quadrangle=0, internalNodes=0)\n\nObsolate, use gmsh script (.geo) instead.\n\nReturn: none\n\nTypes:\n\n``: x\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.getTagForPhysicalName-Tuple{Any}","page":"Functions","title":"LowLevelFEM.getTagForPhysicalName","text":"FEM.getTagForPhysicalName(name)\n\nReturns tags of elements of physical group name.\n\nReturn: tags\n\nTypes:\n\nname: String\ntags: Vector{Integer}\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.initialDisplacement!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.initialDisplacement!","text":"FEM.initialDisplacement!(problem, name, u0; ux=..., uy=..., uz=...)\n\nChanges the displacement values ux, uy and uz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in displacement vector u0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: String \nu0: Vector{Float64}\nux: Float64 \nuy: Float64 \nuz: Float64 \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.initialVelocity!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.initialVelocity!","text":"FEM.initialVelocity!(problem, name, v0; vx=..., vy=..., vz=...)\n\nChanges the velocity values vx, vy and vz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in velocity vector v0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: String \nv0: Vector{Float64}\nvx: Float64 \nvy: Float64 \nvz: Float64 \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.load-Tuple{Any}","page":"Functions","title":"LowLevelFEM.load","text":"FEM.load(name; fx=..., fy=..., fz=...)\n\nGives the intensity of distributed load on name physical group. At least one fx,  fy or fz value have to be given (depending on the dimension of the problem).\n\nReturn: none\n\nTypes:\n\nname: string\nux: double\nuy: double\nuz: double\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.loadVector-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.loadVector","text":"FEM.loadVector(problem, loads)\n\nSolves a load vector of problem. loads is a tuple of name of physical group  name, coordinates fx, fy and fz of the intensity of distributed force. It can solve traction or body force depending on the problem. In case of 2D problems and Line physical group means surface force. In case of 2D problems and Surface physical group means body force. In case of 3D problems and Line physical group means edge force. In case of 3D problems and Surface physical group means surface force. In case of 3D problems and Volume physical group means body force.\n\nReturn: loadVec\n\nTypes:\n\nproblem: Problem\nloads: Vector{Tuple{String, Float64, Float64, Float64}}\nloadVec: Vector\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.massMatrix-Tuple{Any}","page":"Functions","title":"LowLevelFEM.massMatrix","text":"FEM.massMatrix(problem; lumped=...)\n\nSolves the mass matrix of the problem. If lumped is true, solves lumped mass matrix.\n\nReturn: massMat\n\nTypes:\n\nproblem: Problem\nlumped: Boolean\nmassMat: SparseMatrix\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.material-Tuple{Any}","page":"Functions","title":"LowLevelFEM.material","text":"FEM.material(name; E=2.0e5, ν=0.3, ρ=7.85e-9)\n\nReturns a tuple in which name is the name of a physical group,  E is the modulus of elasticity, ν Poisson's ratio and ρ is the mass density.\n\nReturn: mat\n\nTypes:\n\nmat: Tuple(String, Float64, Float64, Float64)\nname: String\nE: double\nν: double\nρ: double\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.nodalAcceleration!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.nodalAcceleration!","text":"FEM.nodalAcceleration!(problem, name, a0; ax=..., ay=..., az=...)\n\nChanges the acceleration values ax, ay and az (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in acceleration vector a0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: String \na0: Vector{Float64}\nax: Float64\nay: Float64\naz: Float64\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.nodalForce!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.nodalForce!","text":"FEM.nodalForce!(problem, name, f0; fx=..., fy=..., fz=...)\n\nChanges the force values fx, fy and fz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in load vector f0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: String \nf0: Vector{Float64}\nfx: Float64 \nfy: Float64 \nfz: Float64 \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.plotOnPath-NTuple{4, Any}","page":"Functions","title":"LowLevelFEM.plotOnPath","text":"FEM.plotOnPath(problem, pathName, field, points; step=..., name=..., visible=...)\n\nLoad a 2D plot on a path into a View in gmsh. field is the number of View in gmsh from which the data of a field is imported. pathName is the name of a physical group which contains a curve. The curve is devided into equal length intervals with number of points points. The field is shown at this points. step is the sequence number of displayed step. If no step is given, shows all the aviable steps as an animation. name is the title of graph and visible is a true or false value to toggle on or off the initial visibility  in gmsh. This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\npathName: String\nfield: Integer\npoints: Integer\nstep: Integer\nname: String\nvisible: Boolean\ntag: Integer\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.showDoFResults-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.showDoFResults","text":"FEM.showDoFResults(problem, q, comp; t=..., name=..., visible=...)\n\nLoads nodal results into a View in gmsh. q is the field to show, comp is the component of the field (\"uvec\", \"ux\", \"uy\", \"uz\", \"vvec\", \"vx\", \"vy\", \"vz\"), t is a vector of time steps (same number of columns as q), name is a title to display and visible is a true or false value to toggle on or off the  initial visibility in gmsh. If q has more columns, then a sequence of results will be shown (eg. as an animation). This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\nq: Vector{Matrix}\ncomp: String\nt: Vector{Float64}\nname: String\nvisible: Boolean\ntag: Integer\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.showStressResults-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.showStressResults","text":"FEM.showStressResults(problem, S, comp; t=..., name=..., visible=..., smooth=...)\n\nLoads stress results into a View in gmsh. S is a stress field to show, comp is the component of the field (\"s\", \"sx\", \"sy\", \"sz\", \"sxy\", \"syz\", \"szx\", \"seqv\"), t is a vector of time steps (same length as the number of stress states), name is a title to display, visible is a true or false value to toggle on or off the initial visibility in gmsh and smooth is a true of false value to toggle smoothing the stress field on or off. If length of t is more than one, then a  sequence of results will be shown (eg. as an animation). This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\nS: StressField\ncomp: String\nt: Vector{Float64}\nname: String\nvisible: Boolean\nsmooth: Boolean\ntag: Integer\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.smallestPeriodTime-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.smallestPeriodTime","text":"FEM.smallestPeriodTime(K, M)\n\nSolves the smallest period of time for a dynamic problem given by stiffness matrix K and the mass matrix M.`\n\nReturn: Δt\n\nTypes:\n\nK: SparseMatrix\nM: SparseMatrix\nΔt: Float64 \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.solveDisplacement-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.solveDisplacement","text":"FEM.solveDisplacement(K, q)\n\nSolves the equation K*q=f for the displacement vector q. K is the stiffness Matrix, q is the load vector.\n\nReturn: q\n\nTypes:\n\nK: SparseMatrix \nf: Vector{Float64} \nq: Vector{Float64}\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.solveStress-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.solveStress","text":"FEM.solveStress(problem, q)\n\nSolves the stress field S from displacement vector q. Stress field is given per elements, so it usually contains jumps at the boundary of elements. Details of mesh is available in problem.\n\nReturn: S\n\nTypes:\n\nproblem: Problem\nq: Vector{Float64}\nS: StressField\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.stiffnessMatrix-Tuple{Any}","page":"Functions","title":"LowLevelFEM.stiffnessMatrix","text":"FEM.stiffnessMatrix(problem)\n\nSolves the stiffness matrix of the problem.\n\nReturn: stiffMat\n\nTypes:\n\nproblem: Problem\nstiffMat: SparseMatrix\n\n\n\n\n\n","category":"method"},{"location":"Functions/#Index","page":"Functions","title":"Index","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"#LowLevelFEM","page":"Introduction","title":"LowLevelFEM","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Solution of a problem in linear elasticity using Finite Element Method consists of solution of the stiffness matrix mathbfK and load vector mathbff, modifying them according to the boundary conditions (getting tildemathbfK and tildemathbff), solving the displacement field mathbfq as the result of system of equations tildemathbfKmathbfq=tildemathbff, solving the stress field from mathbfq and visualize them. The above described steps can be easily performed using the LowLevelFEM package. Each step means a function with the appropriate parameters, while at any step it is possible to perform an arbitrary operation with the quantities calculated in the meantime. For example the strain energy can be solved as U=frac12mathbfq^TmathbfKmathbfq, for which the code is simply U=q'*K*q/2.(see Examples)","category":"page"},{"location":"#Features","page":"Introduction","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Sketching the geometry, making the FE mesh with GMSH.\nSolving problems from linear elasticity,\n2D problems,\nPlane stress,\nPlane strain,\n3D problem (solid body),\nwhich means the creation of the stiffness matrix mathbfK of the problem using arbitrary\nelement types (line, triangle, rectangle, tetrahedron, hexahedron, pyramid, wedge)\napproximation order (up to ten, Lagrange polynomials)\nApplying\ndistributed forces on arbitrary physical groups (see GMSH),\nLines (in 2D: surface force, in 3D: edge force),\nSurfaces (in 2D: body force, in 3D: traction),\nVolumes (in 3D: body force),\nconcentrated forces on nodes, which means the calculation of the load vector mathbff.\nConstraints on physical groups (nodes on points, edges, surfaces and volumes).\nApplying initial conditions on arbitrary points, edges, surfaces, volumes and on combinations of them.\nSolution of static and dynamic (transient with central difference method) problems,\nwhich means the generations of the mass matrix mathbfM.\nDisplaying the results (scalar or vector displacements, scalar or tensor stresses) with GMSH.\nWhen dynamic problems are solved animations are also possible (click on triangleright).\nPlotting arbitrary results on paths.","category":"page"},{"location":"#Planned-features","page":"Introduction","title":"Planned features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"[ ] 2D axisymmetric problem\n[ ] 3D (and 2D) beam structures\n[ ] Shells\n[ ] Giving loads and prescribed displacements with functions\n[ ] MultiPoint Constraint (like MPC184 in Ansys)\n[x] Different material properties on physical groups\n[ ] Contact problems,\n[ ] in 2D,\n[ ] in 3D,\n[ ] with Lagrange multiplier method.\n[ ] Defining and using coordinate systems,\n[ ] cartesian at arbitrary position and arbitrary orientation,\n[ ] cylindrical.\n[ ] Finite deformations.\n[ ] Heat conduction problems,\n[ ] solving conductivity matrix,\n[ ] solving heat capacity matrix.\n[ ] Dynamic transient problems with HHT-α (or Newmark).\n[ ] Linear buckling.\n[ ] Modal analysis (eigenfrequencies, modal shapes).","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Any suggestions are welcome.","category":"page"}]
}
