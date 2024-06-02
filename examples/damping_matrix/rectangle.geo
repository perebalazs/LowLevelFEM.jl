//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("bottom", 7) = {1};
//+
Transfinite Curve {3} = 2 Using Progression 1;
//+
Transfinite Curve {2} = 2 Using Progression 1;
//+
Transfinite Curve {1} = 2 Using Progression 1;
//+
Transfinite Curve {4} = 2 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Mesh 2;