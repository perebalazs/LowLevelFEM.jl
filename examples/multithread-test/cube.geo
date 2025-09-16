//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};

//+
Transfinite Curve {:} = 21 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{1};
//+
Recombine Surface {:};
//+
Mesh.ElementOrder = 2;
//+
Mesh 3;
//+
Physical Surface("top", 13) = {4};
//+
Physical Surface("bottom", 14) = {3};
//+
Physical Surface("left", 15) = {1};
//+
Physical Surface("right", 16) = {2};
//+
Physical Surface("front", 17) = {6};
//+
Physical Surface("rear", 18) = {5};
//+
Physical Volume("body", 19) = {1};
//+
Mesh.SurfaceEdges = 0;
//+
Point(9) = {0.5, 0.5, 0.5, 1.0};
//+
Physical Point("P", 20) = {9};
