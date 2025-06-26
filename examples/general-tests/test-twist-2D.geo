//+
SetFactory("OpenCASCADE");
//+
Disk(1) = {0, 0, 0, 10, 10};
//+
Mesh.ElementOrder=2;
//+
MeshSize {:} = 2;
//+
Mesh 3;
//+
Physical Surface("body", 2) = {1};
//+
Physical Curve("path", 3) = {1};
