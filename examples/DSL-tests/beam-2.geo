//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {200, 0, 0, 1.0};
//+
Point(3) = {200, 100, 0, 1.0};
//+
Point(4) = {300, 100, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Point(5) = {300, 0, 0, 1.0};
//+
Point(6) = {400, 0, 0, 1.0};
//+
Circle(4) = {4, 5, 6};
//+
MeshSize {:} = 10;
Mesh.ElementOrder=1;
Mesh 1;
//+
Physical Curve("beam", 5) = {1, 2, 3, 4};
//+
Physical Point("A", 6) = {1};
//+
Physical Point("B", 7) = {6};
//+
Physical Point("1", 8) = {2};
//+
Physical Point("2", 9) = {3};
//+
Physical Point("3", 10) = {4};
//+
Physical Point("O", 11) = {5};
