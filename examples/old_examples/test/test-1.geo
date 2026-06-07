//+
SetFactory("OpenCASCADE");
Rectangle(1) = {100, 0, 0, 100, 100, 0};
//+
Mesh.ElementOrder=1;
//+
MeshSize {3, 4, 1, 2} = 100;
//+
Transfinite Curve {4, 3, 2, 1} = 2 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Mesh 2;
//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("right", 7) = {2};
//+
Physical Curve("top", 8) = {3};
//+
Physical Curve("bottom", 9) = {1};
//+
Physical Point("leftBottom", 10) = {1};
