//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 10, -2, 10};
//+
Box(2) = {0.0, 0.1, 0.0, 10, 2, 10};
//+
MeshSize {:} = 0.5;
//+
MeshSize {10, 9, 13, 14} = 0.5;
//+
MeshSize {3,4,7,8} = 0.5;
//+
Mesh.ElementOrder=1;
Mesh 3;
//+
Physical Volume("body", 25) = {1, 2};
//+
Physical Surface("bottom", 26) = {3};
//+
Physical Surface("top", 27) = {10};
//+
Physical Surface("master", 28) = {4};
//+
Physical Surface("slave", 29) = {9};
