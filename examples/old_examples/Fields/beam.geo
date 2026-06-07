//+
SetFactory("OpenCASCADE");
Box(1) = {-5, -5, 0, 10, 10, 100};
//+
Transfinite Curve {3, 7, 5, 1} = 41 Using Progression 1;
//+
Transfinite Curve {6, 12, 2, 10, 8, 11, 4, 9} = 5 Using Progression 1;
//+
Transfinite Surface {1:6};
//+
Transfinite Volume{1};
//+
Recombine Surface {1:6};
//+
Mesh.ElementOrder = 1;
//+
Mesh 3;
//+
Physical Volume("body", 13) = {1};
//+
Physical Surface("right", 14) = {6};
//+
Physical Surface("left", 15) = {5};
//+
Physical Surface("top", 16) = {4};
//+
Physical Surface("bottom", 17) = {3};
//+
Physical Surface("front", 18) = {1};
//+
Physical Surface("rear", 19) = {2};
//+
Point(9) = {0, -5, 50, 1.0};
//+
Point(10) = {0, 5, 50, 1.0};
//+
Line(13) = {9, 10};
//+
Physical Curve("path", 20) = {13};
