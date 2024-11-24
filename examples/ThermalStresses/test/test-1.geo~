//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {120, 20, 0, 1.0};
//+
Point(4) = {120, 150, 0, 1.0};
//+
//b = DefineNumber[ 10, Name "Parameters/b" ];
b = 10;
//+
Point(5) = {0, b, 0, 1.0};
//+
Point(6) = {100, b, 0, 1.0};
//+
Point(7) = {120-b, 20, 0, 1.0};
//+
Point(8) = {120-b, 150, 0, 1.0};
//+
Point(9) = {100, 20, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 8};
//+
Line(4) = {8, 7};
//+
Line(5) = {6, 5};
//+
Line(6) = {5, 1};
//+
Circle(7) = {6, 9, 7};
//+
Circle(8) = {2, 9, 3};
//+
Line(9) = {2, 6};
//+
Line(10) = {3, 7};
//+
Curve Loop(1) = {6, 1, 9, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 7, -10, -8};
//+
Plane Surface(2) = {-2};
//+
Curve Loop(3) = {2, 3, 4, -10};
//+
Plane Surface(3) = {3};
//+
Recombine Surface {1, 2, 3};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Curve {7, 8} = 20 Using Progression 1;

Mesh 2;
//+
Physical Surface("body", 11) = {1, 2, 3};
//+
Physical Curve("hob", 12) = {1};
//+
Physical Curve("water", 13) = {5, 7, 4};
//+
Physical Curve("air", 14) = {3, 2, 8};
//+
Physical Point("support", 15) = {1};
