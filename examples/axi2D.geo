//+
Point(1) = {10, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 20, 0, 1.0};
//+
Point(4) = {10, 20, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("supp", 6) = {1};
//+
Physical Curve("load", 7) = {3};

//+
Transfinite Curve {4, 2} = 11 Using Progression 1;
//+
Transfinite Curve {1, 3} = 6 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};


Mesh 2;
