//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 10, 5, 0};
//+
Transfinite Curve {3, 1} = 2 Using Progression 1;
//+
Transfinite Curve {4, 2} = 2 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Recombine Surface {1};
//+
Mesh 2;
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("right", 6) = {2};
//+
Physical Surface("rect", 7) = {1};
//+
Physical Curve("bottom", 8) = {1};
//+
Physical Curve("top", 9) = {3};
