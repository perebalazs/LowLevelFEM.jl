SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
MeshSize {:} = 0.2;
Mesh 3;
Point(9) = {0.5, 0.5, 0.5, 1.0};
Physical Volume("cube", 13) = {1};
Physical Point("P", 14) = {9};
