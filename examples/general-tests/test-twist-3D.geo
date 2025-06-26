//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 100, 0, 0, 10, 2*Pi};
//+
Rotate {{0, 1, 0}, {0, 0, 0}, -Pi/2} {
  Point{1}; Point{2}; Curve{2}; Curve{1}; Curve{3}; Surface{2}; Surface{1}; Surface{3}; Volume{1}; 
}
//+
Physical Volume("body", 7) = {1};
//+
Physical Surface("supp", 8) = {3};
//+
Physical Surface("load", 9) = {2};
//+
Mesh.ElementOrder=2;
//+
MeshSize {1, 2} = 5;
//+
Mesh 3;
