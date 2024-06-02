//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, -1, 0, 100, 2, 0};
//+
Transfinite Curve {4, 2} = 3 Using Progression 1;
//+
Transfinite Curve {3, 1} = 101 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Mesh 2;

//+
Physical Curve("supp", 5) = {2};
//+
Physical Surface("body", 6) = {1};
//+
Point(5) = {0, 0, 0, 1.0};
//+
Point(6) = {100, 0, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Physical Curve("path", 7) = {5};
