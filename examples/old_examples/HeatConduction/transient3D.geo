//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 100, 10, 10};
//+
MeshSize {:} = 3;
//+
Mesh 3;
//+
Point(9) = {0, 5, 5, 1.0};
//+
Point(10) = {100, 5, 5, 1.0};
//+
Line(13) = {9, 10};
//+
Physical Volume("body", 14) = {1};
//+
Physical Surface("left", 15) = {1};
//+
Physical Surface("right", 16) = {2};
//+
Physical Curve("path", 17) = {13};
