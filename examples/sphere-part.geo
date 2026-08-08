// One octant of a sphere for the Stokes theorem tutorial.
// The radius R is set from the Julia notebook.

SetFactory("Built-in");

// Three points on the coordinate axes and the sphere center.
Point(1) = {R, 0, 0, 1.0};
Point(2) = {0, R, 0, 1.0};
Point(3) = {0, 0, R, 1.0};
Point(4) = {0, 0, 0, 1.0};

// The three quarter-circle arcs form the boundary of the spherical patch.
Circle(1) = {1, 4, 2};
Circle(2) = {2, 4, 3};
Circle(3) = {3, 4, 1};
Curve Loop(1) = {2, 3, 1};

// Force the triangular surface onto the sphere centered at Point 4.
Surface(1) = {1} In Sphere {4};

// Connect the three axis points to the center to close the volume.
Line(4) = {4, 1};
Line(5) = {4, 2};
Line(6) = {4, 3};

// Coordinate-plane faces of the spherical octant.
Curve Loop(3) = {5, 2, -6};
Plane Surface(2) = {3};
Curve Loop(4) = {6, 3, -4};
Plane Surface(3) = {4};
Curve Loop(5) = {4, 1, -5};
Plane Surface(4) = {5};

// The spherical patch and the three plane faces bound the volume.
// The plane faces are reversed so every normal points out of the volume.
Surface Loop(1) = {-2, -4, -3, 1};
Volume(1) = {1};

// Use second-order curved tetrahedra and a moderate mesh size.
MeshSize {:} = 1;
Mesh.ElementOrder = 2;
Mesh 3;

// Names used by LowLevelFEM in the tutorial.
Physical Curve("peri", 4) = {1, 2, 3};
Physical Surface("surf", 5) = {1};
Physical Volume("volu", 7) = {1};
Physical Curve("line", 8) = {6, 4, 5};
