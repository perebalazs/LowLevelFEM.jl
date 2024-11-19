//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 100, 10, 0};
//+
Physical Curve("left", 5) = {4};
//+
Physical Curve("right", 6) = {2};
//+
Physical Surface("body", 7) = {1};

Recombine Surface {1};
Transfinite Line {2,4} = 4;
Transfinite Line {1,3} = 31;
Transfinite Surface {1};
Mesh.ElementOrder = 3;

SetName "bending2D";
Mesh 2;
// Mesh.SaveAll=1;
// Save "bending3D.msh";

//+
Point(5) = {0, 5, 0, 1.0};
//+
Point(6) = {100, 5, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Physical Curve("path", 8) = {5};
