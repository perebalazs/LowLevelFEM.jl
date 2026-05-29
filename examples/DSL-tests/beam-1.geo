//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {6000, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Point(3) = {6000, 6000, 0, 1.0};
//+
Line(2) = {2, 3};
//+
MeshSize {:} = 1;
Mesh.ElementOrder=1;
Mesh 1;
//+
Physical Curve("beam", 2) = {1,2};
//+
Physical Point("A", 3) = {1};
//+
Physical Point("B", 4) = {2};
//+
Physical Point("C", 5) = {3};
