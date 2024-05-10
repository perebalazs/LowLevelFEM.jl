
l = 100;
h = 5;
//+
Point(1) = {0, -l, 0, 1.0};
//+
Point(2) = {l, -l, 0, 1.0};
//+
Point(3) = {l, l, 0, 1.0};
//+
Point(4) = {0, l, 0, 1.0};
//+
Point(5) = {l, -h, 0, 1.0};
//+
Point(6) = {l, h, 0, 1.0};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 1};
//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("body", 7) = {1};
//+
Physical Curve("supp", 8) = {5, 6, 1};
//+
Physical Curve("load", 9) = {3};

Mesh 2;
