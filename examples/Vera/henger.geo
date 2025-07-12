// Gmsh project created on Mon Jun 30 10:31:07 2025
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 20, 20, 2*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 20, 40, 2*Pi};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
//+
MeshSize {:} = 2;
//+
Mesh 3;
//+
Physical Volume("body", 7) = {2};
//+
Physical Surface("inner", 8) = {1};
//+
Physical Surface("outer", 9) = {2};
//+
Physical Surface("front", 10) = {3};
//+
Physical Surface("rear", 11) = {4};
