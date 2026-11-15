//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {10, 10, 10, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Mesh 1;
//+
Physical Point("A", 4) = {1};
//+
Physical Point("B", 5) = {2};
//+
Physical Point("C", 6) = {3};
//+
Physical Point("D", 7) = {4};
//+
Physical Curve("beam", 8) = {1, 2, 3};
//+
Physical Curve("AB", 9) = {1};
//+
Physical Curve("BC", 10) = {2};
//+
Physical Curve("CD", 11) = {3};
