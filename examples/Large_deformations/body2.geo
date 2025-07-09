//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 10, 1, 1};
//+
//Coherence;
//+
//MeshSize {:} = 1;
//+
Transfinite Curve {9:12} = 21 Using Progression 1;
Transfinite Curve {1:8} = 5 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{:};
//+
Recombine Surface {:};
//+
Mesh.ElementOrder=1;
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
//+
Point(9) = {5, 0, 0.5, 1.0};
//+
Point(10) = {5, 1, 0.5, 1.0};
//+
Point(11) = {5, 0.5, 0, 1.0};
//+
Point(12) = {5, 0.5, 1, 1.0};
//+
Line(13) = {9, 10};
//+
Line(14) = {12, 11};
//+
Physical Curve("horizontal", 29) = {14};
//+
Physical Curve("vertical", 30) = {13};
//+
Point(13) = {5, 0.3, 0.3, 1.0};
//+
Physical Point("A", 31) = {13};
