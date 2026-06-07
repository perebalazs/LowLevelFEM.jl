//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 100, 10, 0};
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
MeshSize {:} = 2.5;
Mesh.ElementOrder = 2;
Mesh 2;
//+
Point(5) = {0, 5, 0, 1.0};
//+
Point(6) = {100, 5, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Physical Surface("body", 6) = {1};
//+
Physical Curve("left", 7) = {4};
//+
Physical Curve("right", 8) = {2};
//+
Physical Curve("path", 9) = {5};
