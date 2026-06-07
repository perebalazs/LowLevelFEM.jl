//+
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 5, 0, 1, 0, 5, 2*Pi};
//+
Mesh 3;
//+
Physical Volume("wheel", 4) = {1};
