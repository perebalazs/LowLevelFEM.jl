//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
a = 0.;
//+
Point(1) = {-10, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {50, 0, 0, 1.0};
//+
Point(4) = {50, -50, 0, 1.0};
//+
Point(5) = {-50, -50, 0, 1.0};
//+
Point(6) = {-50, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {-5, -6, -1, -2, -3, -4};
//+
Plane Surface(1) = {1};
//+
Point(7) = {-50, 50+a, 0, 1.0};
//+
Point(8) = {-10, 50-Sqrt(50^2-10^2)+a, 0, 1.0};
//+
Point(9) = {10, 50-Sqrt(50^2-10^2)+a, 0, 1.0};
//+
Point(10) = {50, 50+a, 0, 1.0};
//+
Point(11) = {0, 50+a, 0, 1.0};
//+
Circle(7) = {7, 11, 8};
//+
Circle(8) = {8, 11, 9};
//+
Circle(9) = {9, 11, 10};
//+
Line(10) = {7, 11};
//+
Line(11) = {11, 10};
//+
Curve Loop(2) = {7, 8, 9, -11, -10};
//+
Plane Surface(2) = {2};
//+
MeshSize {7, 10, 11, 4, 5, 3, 6} = 10;
//+
MeshSize {1, 2, 8, 9} = 0.2;

Mesh.ElementOrder=2;

Mesh 2;//+
Physical Surface("upper", 12) = {2};
//+
Physical Surface("lower", 13) = {1};
//+
Physical Curve("bottom", 14) = {4};
//+
Physical Curve("left", 15) = {5};
//+
Physical Curve("right", 16) = {3};
//+
Physical Curve("top", 17) = {11, 10};
//+
Physical Curve("master", 18) = {8};
//+
Physical Curve("slave", 19) = {1};
