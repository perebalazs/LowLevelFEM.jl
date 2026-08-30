//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {200, 0, 0, 1.0};
//+
Point(4) = {200, 0, 0, 1.0};
//+
Point(5) = {300, 0, 0, 1.0};
//+
Point(6) = {400, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
MeshSize {:} = 10;
Mesh.ElementOrder=3;
//+
Mesh 1;
//+
Physical Curve("ABC", 5) = {1, 2};
//+
Physical Curve("CDE", 6) = {3, 4};
//+
Physical Point("A", 7) = {1};
//+
Physical Point("B", 8) = {2};
//+
Physical Point("C1", 9) = {3};
//+
Physical Point("C2", 10) = {4};
//+
Physical Point("D", 11) = {5};
//+
Physical Point("E", 12) = {6};
//+
Physical Curve("AB", 13) = {1};
