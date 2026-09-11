//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 10, -2, 10};
//+
Box(2) = {0, 0.1, 0, 10, 2, 10};
//+
MeshSize {12, 11, 10, 9, 14, 13, 15, 16} = 0.7;
//+
MeshSize {4, 3, 2, 1, 8, 7, 6, 5} = 1;
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
