var documenterSearchIndex = {"docs":
[{"location":"index.html#LowLevelFEM.jl-Documentation","page":"Introduction","title":"LowLevelFEM.jl Documentation","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"index.html#Functions","page":"Introduction","title":"Functions","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":" Modules = [LowLevelFEM]","category":"page"},{"location":"index.html#LowLevelFEM.Problem","page":"Introduction","title":"LowLevelFEM.Problem","text":"Problem(; E=..., ν=..., ρ=..., thickness=..., type=...)\n\nA structure containing the most important data of the problem. \n\nname of the model (in gmsh)\ntype of the problem: 3D \"Solid\", \"PlaneStrain\" or \"PlaneStress\"\ndimension of the problem, determined from type\nmaterial constants: Young's modulus, Poisson's ratio, mass density\nthickness of the plate\nnumber of nodes (non)\n\nTypes:\n\nname: string\ntype: string\ndim: Integer\nE: double\nν: double\nρ: double\nthickness: double\nnon: Integer\n\n\n\n\n\n","category":"type"},{"location":"index.html#LowLevelFEM.StressField","page":"Introduction","title":"LowLevelFEM.StressField","text":"StressField(sigma, numElem, nsteps)\n\nA structure containing the data of a stress field. \n\nsigma: vector of ElementNodeData type stress data (see gmsh.jl)\nnumElem: vector of tags of elements\nnsteps: number of stress fields stored in sigma (for animations).\n\nTypes:\n\nsigma: Vector{Matrix{Float64}}\nnumElem: Vector{Integer}\nnsteps: Integer\n\n\n\n\n\n","category":"type"},{"location":"index.html#LowLevelFEM.CDM-NTuple{8, Any}","page":"Introduction","title":"LowLevelFEM.CDM","text":"FEM.CDM(K, M, C, f, u0, v0, T, Δt)\n\nSolves a transient dynamic problem using central difference method (explicit). K is the stiffness Matrix, M is the mass matrix, C is the damping matrix, f is the load vector, u0 is the initial displacement, v0 is the initial velocity, T is the upper bound ot the time intervall (lower bound is zero) and Δt is the time step size. Returns the displacement vectors and velocity vectors in each time step arranged in the columns of the two matrices u and v and a vector t of the time instants used.\n\nReturn: u, v, t\n\nTypes:\n\nK: Matrix\nM: Matrix\nC: Matrix\nf: Vector\nu0: Vector\nv0: Vector\nT: Double \nΔt: Double \nu: Matrix\nv: Matrix\nt: Vector\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.applyBoundaryConditions!-NTuple{4, Any}","page":"Introduction","title":"LowLevelFEM.applyBoundaryConditions!","text":"FEM.applyBoundaryConditions!(problem, stiffMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nstiffMat: Matrix \nloadVec: Vector \nsupports: Tuple(string, double, double, double)\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.applyBoundaryConditions!-NTuple{6, Any}","page":"Introduction","title":"LowLevelFEM.applyBoundaryConditions!","text":"FEM.applyBoundaryConditions!(problem, stiffMat, massMat, dampMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat, mass matrix massMat, damping matrix dampMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nstiffMat: Matrix \nloadVec: Vector \nsupports: Tuple(string, double, double, double)\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.applyBoundaryConditions-NTuple{4, Any}","page":"Introduction","title":"LowLevelFEM.applyBoundaryConditions","text":"FEM.applyBoundaryConditions(problem, stiffMat, loadVec, supports)\n\nApplies displacement boundary conditions supports on a stiffness matrix stiffMat and load vector loadVec. Mesh details are in problem. supports is a tuple of name of physical group and prescribed displacements ux, uy and uz. Creates a new stiffness matrix and load vector.\n\nReturn: stiffMat, loadVec\n\nTypes:\n\nproblem: Problem\nstiffMat: Matrix \nloadVec: Vector \nsupports: Tuple(string, double, double, double)\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.displacementConstraint-Tuple{Any}","page":"Introduction","title":"LowLevelFEM.displacementConstraint","text":"FEM.displacementConstraint(name; ux=..., uy=..., uz=...)\n\nGives the displacement constraints on name physical group. At least one ux,  uy or uz value have to be given (depending on the dimension of the problem).\n\nReturn: none\n\nTypes:\n\nname: string\nux: double\nuy: double\nuz: double\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.generateMesh-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.generateMesh","text":"FEM.generateMesh(...)\n\nGives...\n\nReturn: none\n\nTypes:\n\n``: x\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.getTagForPhysicalName-Tuple{Any}","page":"Introduction","title":"LowLevelFEM.getTagForPhysicalName","text":"FEM.getTagForPhysicalName(name)\n\nReturns tags of elements of physical group name.\n\nReturn: tags\n\nTypes:\n\nname: string\ntags: vector of integers\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.initialDisplacement!-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.initialDisplacement!","text":"FEM.initialDisplacement!(problem, name, u0; ux=..., uy=..., uz=...)\n\nChanges the displacement values ux, uy and uz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in displacement vector u0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \nu0: Vector \nux: Double \nuy: Double \nuz: Double \n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.initialVelocity!-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.initialVelocity!","text":"FEM.initialVelocity!(problem, name, v0; vx=..., vy=..., vz=...)\n\nChanges the velocity values vx, vy and vz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in velocity vector v0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \nv0: Vector \nvx: Double \nvy: Double \nvz: Double \n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.load-Tuple{Any}","page":"Introduction","title":"LowLevelFEM.load","text":"FEM.load(name; fx=..., fy=..., fz=...)\n\nGives the intensity of distributed load on name physical group. At least one fx,  fy or fz value have to be given (depending on the dimension of the problem).\n\nReturn: none\n\nTypes:\n\nname: string\nux: double\nuy: double\nuz: double\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.loadVector-Tuple{Any, Any}","page":"Introduction","title":"LowLevelFEM.loadVector","text":"FEM.loadVector(problem, loads)\n\nSolves a load vector of problem. loads is a tuple of name of physical group  name, coordinates fx, fy and fz of the intensity of distributed force. It can solve traction or body force depending on the problem. In case of 2D problems and Line physical group means surface force. In case of 2D problems and Surface physical group means body force. In case of 3D problems and Line physical group means edge force. In case of 3D problems and Surface physical group means surface force. In case of 3D problems and Volume physical group means body force.\n\nReturn: loadVec\n\nTypes:\n\nproblem: Problem\nloads: Tuple(string, double, double, double)\nloadVec: Vector\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.massMatrix-Tuple{Any}","page":"Introduction","title":"LowLevelFEM.massMatrix","text":"FEM.massMatrix(problem; name, ρ=..., lumped=...)\n\nSolves the mass matrix of the problem. If lumped is true, solves lumped mass matrix.\n\nReturn: massMat\n\nTypes:\n\nproblem: Problem\nname: string - unused\nρ: double - unused\nlumped: boolean\nmassMat: Matrix\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.nodalAcceleration!-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.nodalAcceleration!","text":"FEM.nodalAcceleration!(problem, name, a0; ax=..., ay=..., az=...)\n\nChanges the acceleration values ax, ay and az (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in acceleration vector a0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \na0: Vector \nax: Double \nay: Double \naz: Double \n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.nodalForce!-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.nodalForce!","text":"FEM.nodalForce!(problem, name, f0; fx=..., fy=..., fz=...)\n\nChanges the force values fx, fy and fz (depending on the dimension of the problem) at nodes belonging to physical group name. Original values are in load vector f0.\n\nReturn: none\n\nTypes:\n\nproblem: Problem\nname: string \nf0: Vector \nfx: Double \nfy: Double \nfz: Double \n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.plotOnPath-NTuple{4, Any}","page":"Introduction","title":"LowLevelFEM.plotOnPath","text":"FEM.plotOnPath(problem, pathName, field, points; numOfSteps=..., name=..., visible=...)\n\nLoad a 2D plot on a path into a View in gmsh. field is the number of View in gmsh from which the data of a field is imported. pathName is the name of a physical group which contains a curve. The curve is devided into equal length intervals with number of points points. The field is shown at this points. numOfSteps is the sequence number of steps. name is the title of graph and visible is a true or false value to toggle on or off the initial visibility  in gmsh. This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\npathName: string\nfield: Integer\npoints: Integer\nnumOfStep: Integer\nname: String\nvisible: Boolean\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.showDoFResults-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.showDoFResults","text":"FEM.showDoFResults(problem, q, comp; t=..., name=..., visible=...)\n\nLoads nodal results into a View in gmsh. q is the field to show, comp is the component of the field (\"uvec\", \"ux\", \"uy\", \"uz\", \"vvec\", \"vx\", \"vy\", \"vz\"), t is a vector of time steps (same number of columns as q), name is a title to display and visible is a true or false value to toggle on or off the  initial visibility in gmsh. If q has more columns, then a sequence of results will be shown (eg. as an animation). This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\nq: Vector or Matrix\ncomp: string\nt: Vector\nname: string\nvisible: Boolean\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.showStressResults-Tuple{Any, Any, Any}","page":"Introduction","title":"LowLevelFEM.showStressResults","text":"FEM.showStressResults(problem, S, comp; t=..., name=..., visible=..., smooth=...)\n\nLoads stress results into a View in gmsh. S is a stress field to show, comp is the component of the field (\"s\", \"sx\", \"sy\", \"sz\", \"sxy\", \"syz\", \"szx\"), t is a vector of time steps (same length as the number of stress states), name is a title to display, visible is a true or false value to toggle on or off the initial visibility in gmsh and smooth is a true of false value to toggle smoothing the stress field on or off. If length of t is more than one, then a  sequence of results will be shown (eg. as an animation). This function returns the tag of View.\n\nReturn: tag\n\nTypes:\n\nproblem: Problem\nS: StressField\ncomp: string\nt: Vector\nname: string\nvisible: Boolean\nsmooth: Boolean\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.smallestPeriodTime-Tuple{Any, Any}","page":"Introduction","title":"LowLevelFEM.smallestPeriodTime","text":"FEM.smallestPeriodTime(K, M)\n\nSolves the smallest period of time for a dynamic problem given by stiffness matrix K and the mass matrix M.`\n\nReturn: Δt\n\nTypes:\n\nK: Matrix\nM: Matrix\nΔt: Double \n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.solveDisplacement-Tuple{Any, Any}","page":"Introduction","title":"LowLevelFEM.solveDisplacement","text":"FEM.solveDisplacement(K, q)\n\nSolves the equation K*q=f for the displacement vector q. K is the stiffness Matrix, q is the load vector.\n\nReturn: q\n\nTypes:\n\nK: Matrix \nf: Vector \nq: Vector \n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.solveStress-Tuple{Any, Any}","page":"Introduction","title":"LowLevelFEM.solveStress","text":"FEM.solveStress(problem, q)\n\nSolves the stress field S from displacement vector q. Stress field is given per elements, so it usually contains jumps at the boundary of elements. Details of mesh is available in problem.\n\nReturn: S\n\nTypes:\n\nproblem: Problem\nq: Vector\nS: StressField\n\n\n\n\n\n","category":"method"},{"location":"index.html#LowLevelFEM.stiffnessMatrix-Tuple{Any}","page":"Introduction","title":"LowLevelFEM.stiffnessMatrix","text":"FEM.stiffnessMatrix(problem; name, E=..., ν=...)\n\nSolves the stiffness matrix of the problem.\n\nReturn: stiffMat\n\nTypes:\n\nproblem: Problem\nname: string - unused\nE: double - unused\nν: double - unused\nstiffMat: Matrix\n\n\n\n\n\n","category":"method"},{"location":"index.html#Index","page":"Introduction","title":"Index","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"","category":"page"}]
}
