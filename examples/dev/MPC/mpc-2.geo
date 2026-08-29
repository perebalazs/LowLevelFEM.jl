//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 100, 10, 2*Pi};
//+
Point(3) = {0, 0, 101, 1.0};
//+
MeshSize {:} = 5;
Mesh.ElementOrder=2;
//+
Mesh 3;
//+
Physical Volume("body", 4) = {1};
//+
Physical Surface("left", 5) = {3};
//+
Physical Surface("right", 6) = {2};
//+
Physical Point("remote", 7) = {3};
