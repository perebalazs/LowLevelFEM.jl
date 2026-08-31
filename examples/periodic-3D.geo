//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 10, 10, 100};
//+
Point(9) = {5, 5, 0, 1.0};
//+
Point(10) = {5, 5, 100, 1.0};
//+
MeshSize {:} = 2;
Mesh 3;
//+
Physical Point("P", 13) = {9};
//+
Physical Point("Q", 14) = {10};
//+
Physical Surface("left", 15) = {5};
//+
Physical Surface("right", 16) = {6};
//+
Physical Volume("body", 17) = {1};
