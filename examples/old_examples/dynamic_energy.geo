//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 20, 10, 0};
//+
Transfinite Curve {3, 1} = 161 Using Progression 1;
//+
Transfinite Curve {4, 2} = 81 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

Mesh.ElementOrder = 1;
Mesh 2;

//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("right", 7) = {2};
