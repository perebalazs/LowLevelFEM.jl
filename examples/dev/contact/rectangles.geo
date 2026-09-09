//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 10, -2, 0};
//+
Rectangle(2) = {0, 0.1, 0, 10, 2, 0};
//+
Transfinite Curve {1} = 11 Using Progression 1;
//+
Transfinite Curve {5} = 10 Using Progression 1;
//+
MeshSize {8, 7, 3, 4} = 1;
//+
Mesh.ElementOrder=2;
Mesh 2;
//+
Physical Surface("body", 9) = {2, 1};
//+
Physical Curve("master", 10) = {1};
//+
Physical Curve("slave", 11) = {5};
//+
Physical Curve("bottom", 12) = {3};
//+
Physical Curve("top", 13) = {7};
