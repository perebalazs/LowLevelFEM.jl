//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 100, 100, 100};
//+
Box(2) = {100, 0, 0, 100, 100, 100};
//+
Box(3) = {200, 0, 0, 100, 100, 100};

Coherence;

//+
Physical Volume("mat1", 37) = {1};
//+
Physical Volume("mat2", 38) = {2};
//+
Physical Volume("mat3", 39) = {3};
//+
Physical Surface("supp", 40) = {1};
//+
Physical Surface("load", 41) = {10};

SetName "mat_test";
Mesh 3;