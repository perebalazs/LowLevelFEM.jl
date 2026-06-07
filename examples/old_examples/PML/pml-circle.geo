//+
SetFactory("OpenCASCADE");
//+
e = 1;
//+
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 15, 0, 2*Pi};
//+
Point(3) = {-0.5*e, -0.86/3*e, 0, 1.0};
//+
Point(4) = {0.5*e, -0.86/3*e, 0, 1.0};
//+
Point(5) = {0, 0.86*2/3*e, 0, 1.0};
//+
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 3};
//+
Curve Loop(1) = {4, 5, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1};
//+
Curve Loop(3) = {5, 3, 4};
//+
Plane Surface(2) = {2, 3};
//+
Curve Loop(4) = {2};
//+
Curve Loop(5) = {1};
//+
Plane Surface(3) = {4, 5};
//+
MeshSize {:} = e;
//+
Mesh 2;
//+
Physical Curve("perimeter", 6) = {2};
//+
Physical Surface("surf", 7) = {3, 2, 1};
//+
Physical Surface("source", 8) = {1};
