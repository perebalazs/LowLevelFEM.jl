//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 40, 20, 0};
//+
Physical Surface("body", 5) = {1};
//+
Recombine Surface {1};
//+
Transfinite Curve {3, 4, 2, 1} = 1001 Using Progression 1;
//+
Transfinite Surface {1};
//+
Mesh 2;