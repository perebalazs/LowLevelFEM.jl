//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Mesh 3;

Physical Volume("body", 13) = {1};
