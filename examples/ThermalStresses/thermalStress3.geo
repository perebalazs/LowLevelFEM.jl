//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 100, 5, 5};
//+
MeshSize {:} = 2;
//+
Mesh 3;
//+
Physical Volume("body", 13) = {1};
//+
Physical Surface("left", 14) = {1};
//+
Physical Surface("right", 15) = {2};
//+
Physical Surface("front", 16) = {6};
//+
Physical Surface("bottom", 17) = {3};
