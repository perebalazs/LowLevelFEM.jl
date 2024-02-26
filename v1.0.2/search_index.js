var documenterSearchIndex = {"docs":
[{"location":"Examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"Examples/#2D-Cantilever","page":"Examples","title":"2D Cantilever","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever2D.jl","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using LinearAlgebra, SparseArrays\nimport LowLevelFEM as FEM\nusing LowLevelFEM\n\ngmsh.initialize()\n\ngmsh.open(\"cantilever2D.geo\")\nproblem = FEM.Problem(type=\"PlaneStress\")\n\nsupp = FEM.displacementConstraint(\"supp\", ux=0, uy=0)\nload = FEM.load(\"load\", fy=-1)\n\nK = FEM.stiffnessMatrix(problem)\nf = FEM.loadVector(problem, [load])\n\nFEM.applyBoundaryConditions!(problem, K, f, [supp])\n\nq = FEM.solveDisplacement(K, f)\nS = FEM.solveStress(problem, q)\n\nu = FEM.showDoFResults(problem, q, \"uvec\", name=\"uvec\", visible=false)\nux = FEM.showDoFResults(problem, q, \"ux\", name=\"ux\", visible=false)\nuy = FEM.showDoFResults(problem, q, \"uy\", name=\"uy\", visible=false)\n\ns = FEM.showStressResults(problem, S, \"s\", name=\"σ\", visible=true, smooth=true)\nsx = FEM.showStressResults(problem, S, \"sx\", name=\"σx\", visible=false, smooth=true)\nsy = FEM.showStressResults(problem, S, \"sy\", name=\"σy\", visible=false, smooth=true)\nsxy = FEM.showStressResults(problem, S, \"sxy\", name=\"τxy\", visible=false, smooth=true)\n\nFEM.plotOnPath(problem, \"path\", sx, 100, name=\"σx\", visible=false);\nFEM.plotOnPath(problem, \"path\", sxy, 100, name=\"τxy\", visible=false);\nFEM.plotOnPath(problem, \"path\", ux, 100, name=\"ux\", visible=false);\n\ngmsh.fltk.run()\ngmsh.finalize()","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever2D.geo","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"SetFactory(\"OpenCASCADE\");\n\nRectangle(1) = {0, 0, 0, 100, 10, 0};\n\nPhysical Curve(\"supp\", 5) = {4};\nPhysical Curve(\"load\", 6) = {2};\nPhysical Surface(\"body\", 7) = {1};\n\nRecombine Surface {1};\n\nTransfinite Line {2,4} = 4;\nTransfinite Line {1,3} = 31;\nTransfinite Surface {1};\n\nMesh.ElementOrder = 3;\n\nSetName \"cantilever2D\";\nMesh 2;\n\nPoint(5) = {10, 0, 0, 1.0};\nPoint(6) = {10, 10, 0, 1.0};\nLine(5) = {5, 6};\n\nPhysical Curve(\"path\", 8) = {5};","category":"page"},{"location":"Examples/#3D-Cantilever","page":"Examples","title":"3D Cantilever","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever3D.jl","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using LinearAlgebra, SparseArrays\nimport LowLevelFEM as FEM\nusing LowLevelFEM\n\ngmsh.initialize()\n\ngmsh.open(\"cantilever3D.geo\")\nproblem = FEM.Problem()\n\nsupp = FEM.displacementConstraint(\"supp\", ux=0, uy=0, uz=0)\nload = FEM.load(\"load\", fy=-1)\n\nK = FEM.stiffnessMatrix(problem)\nf = FEM.loadVector(problem, [load])\n\nFEM.applyBoundaryConditions!(problem, K, f, [supp])\n\nq = FEM.solveDisplacement(K, f)\nS = FEM.solveStress(problem, q)\n\nu = FEM.showDoFResults(problem, q, \"uvec\", name=\"uvec\", visible=false)\nux = FEM.showDoFResults(problem, q, \"ux\", name=\"ux\", visible=false)\nuy = FEM.showDoFResults(problem, q, \"uy\", name=\"uy\", visible=false)\nuz = FEM.showDoFResults(problem, q, \"uz\", name=\"uz\", visible=false)\n\ns = FEM.showStressResults(problem, S, \"s\", name=\"σ\", visible=true, smooth=true)\nsx = FEM.showStressResults(problem, S, \"sx\", name=\"σx\", visible=false, smooth=true)\nsy = FEM.showStressResults(problem, S, \"sy\", name=\"σy\", visible=false, smooth=true)\nsz = FEM.showStressResults(problem, S, \"sz\", name=\"σz\", visible=false, smooth=true)\nsxy = FEM.showStressResults(problem, S, \"sxy\", name=\"τxy\", visible=false, smooth=true)\nsyz = FEM.showStressResults(problem, S, \"syz\", name=\"τyz\", visible=false, smooth=true)\nszx = FEM.showStressResults(problem, S, \"szx\", name=\"τzx\", visible=false, smooth=true)\n\nFEM.plotOnPath(problem, \"path\", sx, 100, name=\"σx\", visible=false);\nFEM.plotOnPath(problem, \"path\", sxy, 100, name=\"τxy\", visible=false);\nFEM.plotOnPath(problem, \"path\", ux, 100, name=\"ux\", visible=false);\n\ngmsh.fltk.run()\ngmsh.finalize()","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"cantilever3D.geo","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"SetFactory(\"OpenCASCADE\");\n\nBox(1) = {0, 0, 0, 100, 10, 10};\n\nPhysical Surface(\"supp\", 13) = {1};\nPhysical Surface(\"load\", 14) = {2};\nPhysical Volume(\"body\", 15) = {1};\n\nRecombine Surface {1:6};\n\nTransfinite Line {1:8} = 4;\nTransfinite Line {9:12} = 31;\nTransfinite Surface {1:6};\nTransfinite Volume {1};\n\nMesh.ElementOrder = 3;\n\nSetName \"cantilever3D\";\nMesh 3;\n\nPoint(9) = {10, 0, 5, 1.0};\nPoint(10) = {10, 10, 5, 1.0};\nLine(13) = {9, 10};\n\nPhysical Curve(\"path\", 16) = {13};","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"Functions/#LowLevelFEM.jl","page":"Functions","title":"LowLevelFEM.jl","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Documentation for LowLevelFEM.jl","category":"page"},{"location":"Functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Modules = [LowLevelFEM]","category":"page"},{"location":"Functions/#LowLevelFEM.Problem","page":"Functions","title":"LowLevelFEM.Problem","text":"Problem(; E=..., ν=..., ρ=..., thickness=..., type=...)\n\nA structure containing the most important data of the problem. \n\nname of the model (in gmsh)\ntype of the problem: 3D \"Solid\", \"PlaneStrain\" or \"PlaneStress\"\ndimension of the problem, determined from type\nmaterial constants: Young's modulus, Poisson's ratio, mass density\nthickness of the plate\nnumber of nodes (non)\n\nTypes:\n\nname: string\ntype: string\ndim: Integer\nE: double\nν: double\nρ: double\nthickness: double\nnon: Integer\n\n\n\n\n\n","category":"type"},{"location":"Functions/#LowLevelFEM.StressField","page":"Functions","title":"LowLevelFEM.StressField","text":"StressField(sigma, numElem, nsteps)\n\nA structure containing the data of a stress field. \n\nsigma: vector of ElementNodeData type stress data (see gmsh.jl)\nnumElem: vector of tags of elements\nnsteps: number of stress fields stored in sigma (for animations).\n\nTypes:\n\nsigma: Vector{Matrix{Float64}}\nnumElem: Vector{Integer}\nnsteps: Integer\n\n\n\n\n\n","category":"type"},{"location":"Functions/#LowLevelFEM.CDM-NTuple{8, Any}","page":"Functions","title":"LowLevelFEM.CDM","text":"FEM.CDM(K, M, C, f, u0, v0, T, Δt)\n\nSolves a transient dynamic problem using central difference method (explicit). K is the stiffness Matrix, M is the mass matrix, C is the damping matrix, f is the load vector, u0 is the initial displacement, v0 is the initial velocity, T is the upper bound ot the time intervall (lower bound is zero) and Δt is the time step size. Returns the displacement vectors and velocity vectors in each time step arranged in the columns of the two matrices u and v and a vector t of the time instants used.\n\nReturn: u, v, t\n\nTypes:\n\nK: Matrix\nM: Matrix\nC: Matrix\nf: Vector\nu0: Vector\nv0: Vector\nT: Double \nΔt: Double \nu: Matrix\nv: Matrix\nt: Vector\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.applyBoundaryConditions!-NTuple{4, Any}","page":"Functions","title":"LowLevelFEM.applyBoundaryConditions!","text":"FEM.applyBoundaryConditions!(problem, stiffMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nstiffMat: Matrix \nloadVec: Vector \nsupports: Tuple(string, double, double, double)\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.applyBoundaryConditions!-NTuple{6, Any}","page":"Functions","title":"LowLevelFEM.applyBoundaryConditions!","text":"FEM.applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat, mass matrix massMat, damping matrix dampMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nstiffMat: Matrix \nloadVec: Vector \nsupports: Tuple(string, double, double, double)\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.applyBoundaryConditions-NTuple{4, Any}","page":"Functions","title":"LowLevelFEM.applyBoundaryConditions","text":"FEM.applyBoundaryConditions(problem, stiffMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz. Creates a new stiffness matrix and load vector.\n\nReturn: stiffMat, loadVec\n\nTypes:\n\nproblem: Problem\nstiffMat: Matrix \nloadVec: Vector \nsupports: Tuple(string, double, double, double)\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.displacementConstraint-Tuple{Any}","page":"Functions","title":"LowLevelFEM.displacementConstraint","text":"FEM.displacementConstraint(name; ux=..., uy=..., uz=...)\n\nGives the displacement constraints on name physical group. At least one ux,  uy or uz value have to be given (depending on the dimension of the problem).\n\nReturn: none\n\nTypes:\n\nname: string\nux: double\nuy: double\nuz: double\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.generateMesh-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.generateMesh","text":"FEM.generateMesh(...)\n\nGives...\n\nReturn: none\n\nTypes:\n\n``: x\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.getTagForPhysicalName-Tuple{Any}","page":"Functions","title":"LowLevelFEM.getTagForPhysicalName","text":"FEM.getTagForPhysicalName(name)\n\nReturns tags of elements of physical group name.\n\nReturn: tags\n\nTypes:\n\nname: string\ntags: vector of integers\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.initialDisplacement!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.initialDisplacement!","text":"FEM.initialDisplacement!(problem, name, u0; ux=..., uy=..., uz=...)\n\nChanges the displacement values ux, uy and uz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in displacement vector u0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \nu0: Vector \nux: Double \nuy: Double \nuz: Double \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.initialVelocity!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.initialVelocity!","text":"FEM.initialVelocity!(problem, name, v0; vx=..., vy=..., vz=...)\n\nChanges the velocity values vx, vy and vz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in velocity vector v0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \nv0: Vector \nvx: Double \nvy: Double \nvz: Double \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.load-Tuple{Any}","page":"Functions","title":"LowLevelFEM.load","text":"FEM.load(name; fx=..., fy=..., fz=...)\n\nGives the intensity of distributed load on name physical group. At least one fx,  fy or fz value have to be given (depending on the dimension of the problem).\n\nReturn: none\n\nTypes:\n\nname: string\nux: double\nuy: double\nuz: double\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.loadVector-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.loadVector","text":"FEM.loadVector(problem, loads)\n\nSolves a load vector of problem. loads is a tuple of name of physical group  name, coordinates fx, fy and fz of the intensity of distributed force. It can solve traction or body force depending on the problem. In case of 2D problems and Line physical group means surface force. In case of 2D problems and Surface physical group means body force. In case of 3D problems and Line physical group means edge force. In case of 3D problems and Surface physical group means surface force. In case of 3D problems and Volume physical group means body force.\n\nReturn: loadVec\n\nTypes:\n\nproblem: Problem\nloads: Tuple(string, double, double, double)\nloadVec: Vector\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.massMatrix-Tuple{Any}","page":"Functions","title":"LowLevelFEM.massMatrix","text":"FEM.massMatrix(problem; name, ρ=..., lumped=...)\n\nSolves the mass matrix of the problem. If lumped is true, solves lumped mass matrix.\n\nReturn: massMat\n\nTypes:\n\nproblem: Problem\nname: string - unused\nρ: double - unused\nlumped: boolean\nmassMat: Matrix\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.nodalAcceleration!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.nodalAcceleration!","text":"FEM.nodalAcceleration!(problem, name, a0; ax=..., ay=..., az=...)\n\nChanges the acceleration values ax, ay and az (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in acceleration vector a0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \na0: Vector \nax: Double \nay: Double \naz: Double \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.nodalForce!-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.nodalForce!","text":"FEM.nodalForce!(problem, name, f0; fx=..., fy=..., fz=...)\n\nChanges the force values fx, fy and fz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in load vector f0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \nf0: Vector \nfx: Double \nfy: Double \nfz: Double \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.plotOnPath-NTuple{4, Any}","page":"Functions","title":"LowLevelFEM.plotOnPath","text":"FEM.plotOnPath(problem, pathName, field, points; numOfSteps=..., name=..., visible=...)\n\nLoad a 2D plot on a path into a View in gmsh. field is the number of View in gmsh from which the data of a field is imported. pathName is the name of a physical group which contains a curve. The curve is devided into equal length intervals with number of points points. The field is shown at this points. numOfSteps is the sequence number of steps. name is the title of graph and visible is a true or false value to toggle on or off the initial visibility  in gmsh. This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\npathName: string\nfield: Integer\npoints: Integer\nnumOfStep: Integer\nname: String\nvisible: Boolean\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.showDoFResults-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.showDoFResults","text":"FEM.showDoFResults(problem, q, comp; t=..., name=..., visible=...)\n\nLoads nodal results into a View in gmsh. q is the field to show, comp is the component of the field (\"uvec\", \"ux\", \"uy\", \"uz\", \"vvec\", \"vx\", \"vy\", \"vz\"), t is a vector of time steps (same number of columns as q), name is a title to display and visible is a true or false value to toggle on or off the  initial visibility in gmsh. If q has more columns, then a sequence of results will be shown (eg. as an animation). This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\nq: Vector or Matrix\ncomp: string\nt: Vector\nname: string\nvisible: Boolean\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.showStressResults-Tuple{Any, Any, Any}","page":"Functions","title":"LowLevelFEM.showStressResults","text":"FEM.showStressResults(problem, S, comp; t=..., name=..., visible=..., smooth=...)\n\nLoads stress results into a View in gmsh. S is a stress field to show, comp is the component of the field (\"s\", \"sx\", \"sy\", \"sz\", \"sxy\", \"syz\", \"szx\"), t is a vector of time steps (same length as the number of stress states), name is a title to display, visible is a true or false value to toggle on or off the initial visibility in gmsh and smooth is a true of false value to toggle smoothing the stress field on or off. If length of t is more than one, then a  sequence of results will be shown (eg. as an animation). This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\nS: StressField\ncomp: string\nt: Vector\nname: string\nvisible: Boolean\nsmooth: Boolean\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.smallestPeriodTime-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.smallestPeriodTime","text":"FEM.smallestPeriodTime(K, M)\n\nSolves the smallest period of time for a dynamic problem given by stiffness matrix K and the mass matrix M.`\n\nReturn: Δt\n\nTypes:\n\nK: Matrix\nM: Matrix\nΔt: Double \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.solveDisplacement-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.solveDisplacement","text":"FEM.solveDisplacement(K, q)\n\nSolves the equation K*q=f for the displacement vector q. K is the stiffness Matrix, q is the load vector.\n\nReturn: q\n\nTypes:\n\nK: Matrix \nf: Vector \nq: Vector \n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.solveStress-Tuple{Any, Any}","page":"Functions","title":"LowLevelFEM.solveStress","text":"FEM.solveStress(problem, q)\n\nSolves the stress field S from displacement vector q. Stress field is given per elements, so it usually contains jumps at the boundary of elements. Details of mesh is available in problem.\n\nReturn: S\n\nTypes:\n\nproblem: Problem\nq: Vector\nS: StressField\n\n\n\n\n\n","category":"method"},{"location":"Functions/#LowLevelFEM.stiffnessMatrix-Tuple{Any}","page":"Functions","title":"LowLevelFEM.stiffnessMatrix","text":"FEM.stiffnessMatrix(problem; name, E=..., ν=...)\n\nSolves the stiffness matrix of the problem.\n\nReturn: stiffMat\n\nTypes:\n\nproblem: Problem\nname: string - unused\nE: double - unused\nν: double - unused\nstiffMat: Matrix\n\n\n\n\n\n","category":"method"},{"location":"Functions/#Index","page":"Functions","title":"Index","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"#LowLevelFEM","page":"Introduction","title":"LowLevelFEM","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Solution of a problem in linear elasticity using Finite Element Method consists of solution of the stiffness matrix mathbfK and load vector mathbff, modifying them according to the boundary conditions (getting tildemathbfK and tildemathbff), solving the displacement field mathbfq as the result of system of equations tildemathbfKmathbfq=tildemathbff, solving the stress field from mathbfq and visualize them. The above described steps can be easily performed using the LowLevelFEM package. Each step means a function with the appropriate parameters, while at any step it is possible to perform an arbitrary operation with the quantities calculated in the meantime. For example the strain energy can be solved as U=frac12mathbfq^TmathbfKmathbfq, for which the code is simply U=q'*K*q/2.(see Examples)","category":"page"},{"location":"#Features","page":"Introduction","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Sketching the geometry, making the FE mesh with GMSH.\nSolving problems from linear elasticity,\n2D problems,\nPlane stress,\nPlane strain,\n3D problem (solid body),\nwhich means the creation of the stiffness matrix mathbfK of the problem using arbitrary\nelement types (line, triangle, rectangle, tetrahedron, hexahedron, pyramid, wedge)\napproximation order (up to ten, Lagrange polynomials)\nApplying\ndistributed forces on arbitrary physical groups (see GMSH),\nLines (in 2D: surface force, in 3D: edge force),\nSurfaces (in 2D: body force, in 3D: traction),\nVolumes (in 3D: body force),\nconcentrated forces on nodes, which means the calculation of the load vector mathbff.\nConstraints on physical groups (nodes on points, edges, surfaces and volumes).\nApplying initial conditions on arbitrary points, edges, surfaces, volumes and on combinations of them.\nSolution of static and dynamic (transient with central difference method) problems,\nwhich means the generations of the mass matrix mathbfM.\nDisplaying the results (scalar or vector displacements, scalar or tensor stresses) with GMSH.\nWhen dynamic problems are solved animations are also possible (click on triangleright).\nPlotting arbitrary results on paths.","category":"page"},{"location":"#Planned-features","page":"Introduction","title":"Planned features","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"[ ] 2D axisymmetric problem\n[ ] 3D (and 2D) beam structures\n[ ] Shells\n[ ] Giving loads and prescribed displacements with functions\n[ ] Different material properties on physical groups\n[ ] Contact problems,\n[ ] in 2D,\n[ ] in 3D,\n[ ] with Lagrange multiplier method.\n[ ] Using different materials within a model.\n[ ] Defining and using coordinate systems,\n[ ] cartesian at arbitrary position and arbitrary orientation,\n[ ] cylindrical.\n[ ] Finite deformations.\n[ ] Heat conduction problems,\n[ ] solving conductivity matrix,\n[ ] solving heat capacity matrix.\n[ ] Dynamic transient problems with HHT-α (or Newmark).\n[ ] Linear buckling.\n[ ] Modal analysis (eigenfrequencies, modal shapes).","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Any suggestions are welcome.","category":"page"}]
}