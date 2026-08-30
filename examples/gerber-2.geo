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
Point(7) = {400, 0, 0, 1.0};
//+
Point(8) = {500, 0, 0, 1.0};
//+
Point(9) = {600, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 9};
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
Physical Curve("EFG", 7) = {5, 6};
//+
Physical Point("A", 8) = {1};
//+
Physical Point("B", 9) = {2};
//+
Physical Point("C1", 10) = {3};
//+
Physical Point("C2", 11) = {4};
//+
Physical Point("D", 12) = {5};
//+
Physical Point("E1", 13) = {6};
//+
Physical Point("E2", 14) = {7};
//+
Physical Point("F", 15) = {8};
//+
Physical Point("G", 16) = {9};
//+
Physical Curve("AB", 17) = {1};
