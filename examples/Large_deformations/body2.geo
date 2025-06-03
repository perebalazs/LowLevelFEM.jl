//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Coherence;
//+
//MeshSize {:} = 1;
//+
Transfinite Curve {:} = 3 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{:};
//+
Recombine Surface {:};
//+
Mesh 3;
//+
Physical Volume("body", 21) = {1};
//+
Physical Surface("left", 23) = {1};
//+
Physical Surface("right", 24) = {2};
//+
Physical Surface("bottom", 25) = {3};
//+
Physical Surface("top", 26) = {4};
//+
Physical Surface("front", 27) = {6};
//+
Physical Surface("rear", 28) = {5};
