//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Point(5) = {1.5, 0.5, 0, 1.0};
//+
Transfinite Curve {3, 4, 1, 2} = 2 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Recombine Surface {1};
//+
Mesh 2;
//+
Physical Surface("body", 5) = {1};
//+
Physical Point("remote", 6) = {5};
//+
Physical Curve("right", 7) = {2};//+
Physical Curve("left", 8) = {4};
