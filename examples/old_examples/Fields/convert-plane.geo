//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Rectangle(2) = {1, 0, 0, 1, 1, 0};
//+
Coherence;
//+
Transfinite Curve {1, 2, 3, 4, 5, 6, 7} = 3 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Recombine Surface {1, 2};
//+
Mesh 2;
//+
Physical Surface("plane1", 8) = {1};
//+
Physical Surface("plane2", 9) = {2};
