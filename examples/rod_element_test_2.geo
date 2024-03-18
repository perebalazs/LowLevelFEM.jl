b = 1e-1;

//+
Point(1) = {10, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, b, 0, 1.0};
//+
Point(4) = {10, b, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {3, 1, 4, 2} = 2 Using Progression 1;
//+
Transfinite Surface {1};
//+
Physical Surface("body", 5) = {1};
//+
Recombine Surface {1};

Mesh 2;
