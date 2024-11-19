//+
SetFactory("OpenCASCADE");
Rectangle(1) = {100, 0, 0, 200, 100, 0};
//+
Physical Surface("wall", 5) = {1};
//+
Physical Curve("inner", 6) = {4};
//+
Physical Curve("outer", 7) = {2};
//+
Mesh.ElementOrder = 5;
//+
Mesh 2;
//+
Point(5) = {100, 50, 0, 1.0};
//+
Point(6) = {300, 50, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Physical Curve("path", 8) = {5};
