//+
Point(1) = {0, 0, 0, 15.0};
//+
Point(2) = {100, 0, 0, 15.0};
//+
Point(3) = {100, 50, 0, 15.0};
//+
Point(4) = {50, 50, 0, 0.5};
//+
Point(5) = {50, 100, 0, 15.0};
//+
Point(6) = {0, 100, 0, 15.0};
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
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("fix", 7) = {5};
//+
Physical Curve("load", 8) = {2};
//+
SetName "Lshape";
Mesh.ElementOrder = 4;
Mesh.HighOrderOptimize = 1;
Mesh 2;
//+
Point(7) = {0, 0, 0, 1.0};
//+
Point(8) = {50, 50, 0, 1.0};
//+
Line(7) = {7, 8};
//+
Physical Curve("path", 9) = {7};
//+
//Mesh.SaveAll=1;
//Save "LshapedPlate.msh";
