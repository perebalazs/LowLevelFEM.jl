// Two square plates sharing curve 2.
SetFactory("Built-in");
n = 1;

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Point(5) = {1, 0, 1};
Point(6) = {1, 1, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, -2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Surface("horizontal") = {1};
Physical Surface("vertical") = {2};
Physical Curve("joint") = {2};
Physical Curve("clamp") = {4};
Physical Curve("load") = {6};
Physical Surface("plates") = {1,2};
Physical Point("A") = {5};
Physical Point("B") = {6};

Transfinite Curve {1:7} = n + 1;
Transfinite Surface {1, 2};
Recombine Surface {1, 2};
Mesh.ElementOrder = 1;
Mesh 2;
//+
//+
