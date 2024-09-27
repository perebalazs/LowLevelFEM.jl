//+
Point(1) = {10, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 20, 0, 1.0};
//+
Point(4) = {10, 20, 0, 1.0};
//+
Point(5) = {10, 10, 0, 1.0};
//+
Point(6) = {20, 10, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 6};
//+
Line(3) = {6, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 1};
//+
Line(7) = {5, 6};
//+
Curve Loop(1) = {1, 2, -7, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 3, 4, 5};
//+
Plane Surface(2) = {2};
//+
Physical Curve("supp", 8) = {1};
//+
Physical Curve("load", 9) = {4};
//+
Physical Surface("body1", 10) = {1};
//+
Physical Surface("body2", 11) = {2};
//+
//Recombine Surface {1, 2};
//+
Transfinite Curve {6, 5, 3, 2, 1, 4, 7} = 3 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};

Mesh 2;
