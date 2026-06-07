//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 200, 100, 0};
//+
Point(5) = {100, 0, 0, 1.0};
//+
Point(6) = {100, 100, 0, 1.0};
//+
Line(5) = {5, 6};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{5}; Delete; }
//+
Physical Surface("Fe", 12) = {1};
//+
Physical Surface("Al", 13) = {2};
//+
Physical Curve("supp", 14) = {6};
//+
Physical Curve("load", 15) = {10};
//+
//MeshSize {1, 6, 5, 4, 3, 2} = 10;
MeshSize {:} = 10;


Mesh 2;
