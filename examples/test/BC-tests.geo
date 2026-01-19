//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 10, 10, 10};
//+
Box(2) = {10, 0, 0, 10, 10, 10};
//+
Coherence;
//+
Recombine Surface {:};
//+
Transfinite Curve {:} = 2 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume {:};
//+
Mesh 3;
//+
Physical Surface("rear", 21) = {5};
//+
Physical Surface("left", 22) = {1};
//+
Physical Surface("bottom", 23) = {3};
//+
Physical Surface("front", 24) = {6};
//+
Physical Surface("right", 25) = {2};
//+
Physical Surface("top", 26) = {4};
//+
Physical Volume("volu", 27) = {1};
