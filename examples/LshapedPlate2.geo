R=1;
es = 15;
//+
Point(1) = {0, 0, 0, es};
//+
Point(2) = {100, 0, 0, es};
//+
Point(3) = {100, 50, 0, es};
//+
Point(4) = {50+R, 50, 0, R*(es/15)};
//+
Point(5) = {50, 50+R, 0, R*(es/15)};
//+
Point(6) = {50, 100, 0, es};
//+
Point(7) = {0, 100, 0, es};
//+
Point(8) = {50+R, 50+R, 0, 0.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Circle(4) = {4, 8, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};
//+
Physical Curve("fix", 8) = {6};
//+
Physical Curve("load", 9) = {2};
//+
Physical Surface("body", 11) = {1};
//+
SetName "Lshape";
Mesh.ElementOrder = 1;
Mesh.HighOrderOptimize = 1;
Mesh 2;
//+
Point(9) = {0, 0, 0, 1.0};
//+
Point(10) = {50+0.415*R, 50+0.415*R, 0, 1.0};
//+
Line(8) = {9, 10};
//+
Physical Curve("path", 10) = {8};
//+
//Mesh.SaveAll=1;
//Save "LshapedPlate.msh";





