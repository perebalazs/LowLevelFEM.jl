//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 10, 10, 10};
//+
Box(2) = {10, 0, 0, 10, 10, 10};
//+
Coherence;
//+
Transfinite Curve {:} = 2 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{:};
//+
Recombine Surface {:};
//+
Mesh 3;
//+
Physical Volume("body1", 21) = {1};
//+
Physical Volume("body2", 22) = {2};
