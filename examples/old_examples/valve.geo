//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {5, 0, 0, 1.0};
//+
Point(3) = {5, 5, 0, 1.0};
//+
Point(4) = {5, 100, 0, 1.0};
//+
Point(5) = {0, 100, 0, 1.0};
//+
Point(6) = {20, 0, 0, 1.0};
//+
Point(7) = {20, 5, 0, 1.0};
//+
Line(1) = {5, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {2, 6};
//+
Line(7) = {6, 7};
//+
Point(8) = {15, 5, 0, 1.0};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 3};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 9, -3, 6, 7};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {1} = 101 Using Progression 1;
//+
Transfinite Curve {4} = 96 Using Progression 1;
//+
Transfinite Curve {2, 5, 3, 7, 8} = 6 Using Progression 1;
//+
Transfinite Curve {6} = 16 Using Progression 1;
//+
Transfinite Curve {9} = 11 Using Progression 1;
//+
Transfinite Surface {2} = {2, 6, 7, 3};
//+
Transfinite Surface {1} = {1, 2, 4, 5};
//+
Recombine Surface {1, 2};
//+
Mesh 2;

//+
Physical Surface("body", 10) = {1, 2};
//+
Physical Curve("supp", 11) = {8};
