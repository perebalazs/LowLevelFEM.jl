//+
SetFactory("OpenCASCADE");
Rectangle(1) = {1, 0, 0, 10, 10, 0};


//+
Transfinite Curve {4, 2} = 1000 Using Progression 1;
//+
Transfinite Curve {3, 1} = 1000 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};

Mesh.ElementOrder = 1;
Mesh 2;
//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("supp", 6) = {1};
//+
Physical Curve("load", 7) = {3};
