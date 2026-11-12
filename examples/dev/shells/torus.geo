// Torus geometry for the Gauss divergence theorem tutorial.
// The major radius R and minor radius r are set from the Julia notebook.

SetFactory("OpenCASCADE");

// Create a complete solid torus centered at the origin.
Torus(1) = {0, 0, 0, R, r, 2*Pi};

// Use third-order curved tetrahedra to represent the curved boundary accurately.
Mesh.ElementOrder = 3;
MeshSize {:} = 1;
Mesh 3;

// Names used by LowLevelFEM for the volume and its closed boundary.
Physical Volume("volu", 3) = {1};
Physical Surface("surf", 4) = {1};
