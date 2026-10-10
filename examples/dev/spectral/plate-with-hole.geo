//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 10};
//+
Point(2) = {200, 0, 0, 10};
//+
Point(3) = {200, 100, 0, 10};
//+
Point(4) = {0, 100, 0, 10};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {100, 50, 0, 20, 0, 2*Pi};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
MeshSize {5} = 1;
Mesh.ElementOrder=4;
//+
Physical Surface("body", 6) = {1};
//+
Physical Curve("left", 7) = {4};
//+
Physical Curve("right", 8) = {2};
//+
Recombine Surface {1};
//+
Mesh 2;
