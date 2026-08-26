//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {5, 0, 0, 1.0};
//+
Point(3) = {5, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {2.5, 0.5, 0, 0.3, 0, 2*Pi};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
MeshSize {4, 1, 3, 2} = 0.1;
MeshSize {5} = 0.01;

Mesh.ElementOrder=2;
Mesh 2;

//+
Physical Surface("pipe", 6) = {1};
//+
Physical Curve("left", 7) = {4};
//+
Physical Curve("right", 8) = {2};
//+
Physical Curve("bottom", 9) = {1};
//+
Physical Curve("top", 10) = {3};
//+
Physical Point("reference", 11) = {2};
//+
Physical Curve("obstacle", 12) = {5};
