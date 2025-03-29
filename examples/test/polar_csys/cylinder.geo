//+
SetFactory("OpenCASCADE");
Cylinder(1) = {0, 0, 0, 0, 0, 100, 10, 2*Pi};
//+
MeshSize {:} = 5;
//+
Mesh.ElementOrder=3;
//+
Mesh 3;
//+
Physical Volume("body", 4) = {1};
//+
Physical Surface("left", 5) = {3};
//+
Physical Surface("right", 6) = {2};
//+
Point(13) = {0, 0, 0.1, 1.0};
//+
Point(14) = {10, 0, 0.1, 1.0};
//+
Point(15) = {5*Sqrt(2), 5*Sqrt(2), 0.1, 1.0};
//+
Line(14) = {13, 14};
//+
Line(15) = {13, 15};
//+
Physical Curve("path1", 7) = {14};
//+
Physical Curve("path2", 8) = {15};
//+
Physical Point("A", 9) = {2};
