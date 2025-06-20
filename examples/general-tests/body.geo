//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 50, 1, 1};
//+
Box(2) = {50, 0, 0, 50, 1, 1};
//+
Coherence;
//+
MeshSize {:} = 0.3;
//+
//Transfinite Curve {:} = 2 Using Progression 1;
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
//+
Physical Surface("left", 23) = {1};
//+
Physical Surface("right", 24) = {7};
