//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 2, 1, 1};
//+
//Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} {
//  Point{7}; Point{5}; Point{3}; Point{1}; Point{8}; Point{6}; Point{4}; Point{2}; Curve{6}; Curve{12}; Curve{2}; Curve{10}; Curve{7}; Curve{8}; Curve{5}; Curve{3}; Curve{11}; Curve{9}; Curve{4}; Curve{1}; Surface{6}; Surface{2}; Surface{4}; Surface{5}; Surface{3}; Surface{1}; Volume{1}; 
//}
//+
//Coherence;
//+
//MeshSize {:} = 10;
//+
Transfinite Curve {:} = 4 Using Progression 1;
Transfinite Curve {9:12} = 3 Using Progression 1;
//+
Transfinite Surface {:};
//+
Transfinite Volume{:};
//+
Mesh.ElementOrder=2;
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
