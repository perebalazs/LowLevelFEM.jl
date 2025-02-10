//+
a=10;
l=100;
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, l, a/2, a};
Box(2) = {0, a/2, 0, l, a/2, a};
//+
Coherence;
//+
Transfinite Curve {10, 12, 9, 20, 11, 19} = 31 Using Progression 1;
//+
Transfinite Curve {1, 3, 14, 17, 7, 5} = 5 Using Progression 1;
//+
Transfinite Curve {4, 15, 2, 13, 18, 8, 6, 16} = 3 Using Progression 1;
//+
Transfinite Surface {7, 1, 11, 6, 9, 4, 3, 8, 10, 2, 5};
//+
Transfinite Volume{1};
//+
Transfinite Volume{2};
//+
Recombine Surface {7, 1, 11, 6, 9, 4, 3, 8, 10, 2, 5};

Mesh.ElementOrder = 2;
//Mesh.SecondOrderIncomplete = 1;
//+
Mesh 3;
//+
Physical Curve("supp", 21) = {3};
//+
Physical Curve("load", 22) = {7};
//+
Physical Volume("body", 23) = {1, 2};
//+
Physical Surface("left", 24) = {1, 7};
//+
Physical Surface("right", 25) = {2, 8};
