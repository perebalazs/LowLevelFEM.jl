//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {2, 1, 0, 1.0};
//+
Point(4) = {0, 1, 1, 1.0};
//+
Line(1) = {1, 3};
Line(2) = {2, 3};
Line(3) = {3, 4};
//+
Transfinite Curve {1,2,3} = 2 Using Progression 1;
//+
Mesh 1;
Physical Curve("rod", 2) = {1,2,3};
//+
Physical Point("1", 4) = {1};
//+
Physical Point("2", 5) = {2};
//+
Physical Point("3", 6) = {3};
//+
Physical Point("4", 7) = {4};
