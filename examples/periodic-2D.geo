//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 10, 1, 0};
//+
Point(5) = {0, 0.5, 0, 1.0};
//+
Point(6) = {10, 0.5, 0, 1.0};
//+
MeshSize {:} = 0.2;
//+
Periodic Curve {4}={-2};
//+
Mesh 2;
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("right", 6) = {2};
//+
Physical Surface("body", 7) = {1};
//+
Physical Point("P", 8) = {5};
//+
Physical Point("Q", 9) = {6};