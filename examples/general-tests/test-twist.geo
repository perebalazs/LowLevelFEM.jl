//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 100, 0, 0, 10, 2*Pi};
//+
Physical Volume("body", 7) = {1};
//+
Physical Surface("supp", 8) = {3};
//+
Physical Surface("load", 9) = {2};
//+
Mesh.ElementOrder=1;
//+
MeshSize {1, 2} = 5;
//+
Mesh 3;
