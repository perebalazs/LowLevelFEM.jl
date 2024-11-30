//+
SetFactory("OpenCASCADE");
//+
w=10;
c=250;
f=50;
//+
Rectangle(1) = {-300, -110, 0, 300, 2610, 0};
//+
Rectangle(2) = {4000, -110, 0, 300, 2610, 0};
//+
Rectangle(3) = {-300, 2500, 0, 4600, 200, 0};
//+
Rectangle(4) = {-300, 2700, 0, 4600, c, 0};
//+
Rectangle(5) = {0, -110, 0, 4000, f, 0};
//+
Rectangle(6) = {0, -110+f, 0, 4000, 60, 0};
//+
Rectangle(7) = {-500, -410, 0, 5000, 300, 0};
//+
Rectangle(8) = {-300-w, -110, 0, w, 2810+c, 0};
//+
Rectangle(9) = {4300, -110, 0, w, 2810+c, 0};
//+
Coherence;

//+
MeshSize {:} = 50;

//+
MeshSize {7, 12, 11, 2, 3, 4} = 17;


//+
Recombine Surface {:};
Mesh.ElementOrder=2;

Mesh 2;
//+
Physical Surface("beton", 34) = {7, 3, 6};
//+
Physical Surface("tegla", 35) = {2, 1};
//+
Physical Surface("hungarocell", 36) = {9, 8, 5};
//+
Physical Surface("kgyapot", 37) = {4};
//+
Physical Curve("belsofalb", 38) = {4};
//+
Physical Curve("belsofalj", 43) = {10};
//+
Physical Curve("padlo", 39) = {22};
//+
Physical Curve("plafon", 40) = {13};
//+
Physical Curve("levego", 41) = {32, 31, 18, 30, 33, 25, 28};
//+
Physical Curve("talaj", 42) = {24, 23, 29};

