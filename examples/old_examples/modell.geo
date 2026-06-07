// Gmsh project created on Sat Feb 25 10:36:46 2023

lc = 2;

Point(1) = {0, 0, 0, lc};
Point(2) = {2, 0, 0, lc};
Point(3) = {2, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {1, 1, 0, lc};
Point(6) = {1, 0, 0, lc};
Point(7) = {1, 0.5, 0, lc};

Line(1) = {1,6};
Line(2) = {2,6};
Line(3) = {5,4};
Line(4) = {4,1};
Line(5) = {5,3};
Line(6) = {3,2};
Line(7) = {6,5};

Curve Loop(1) = {1,7,3,4};
Curve Loop(2) = {7,2,6,5};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Recombine Surface {1};


Mesh 2;//+
Physical Surface("body", 8) = {1, 2};
//+
Physical Curve("supp", 9) = {4};
//+
Physical Curve("load", 10) = {6};
