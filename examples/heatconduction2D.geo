//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 100, 10, 0};
//+
Physical Curve("supp", 5) = {4};
//+
Physical Curve("load", 6) = {2};
//+
Physical Surface("body", 7) = {1};

Recombine Surface {1};
Transfinite Line {2,4} = 4;
Transfinite Line {1,3} = 31;
Transfinite Surface {1};
Mesh.ElementOrder = 1;

SetName "bending2D";
Mesh 2;
// Mesh.SaveAll=1;
// Save "bending3D.msh";

//+
Point(5) = {10, 0, 0, 1.0};
//+
Point(6) = {10, 10, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Physical Curve("path", 8) = {5};
