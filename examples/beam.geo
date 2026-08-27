//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 50, 0, 1.0};
//+
Point(4) = {150, 50, 0, 1.0};
//+
Point(5) = {200, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
n = DefineNumber[ 1, Name "Parameters/n" ];
//+
p = DefineNumber[ 4, Name "Parameters/p" ];
//+
Transfinite Curve {1} = 2*n+1 Using Progression 1;
//+
Transfinite Curve {2, 3, 4} = n+1 Using Progression 1;
//+
Mesh.ElementOrder=p;
//+
Mesh 1;
//+
Physical Curve("beam", 5) = {1, 2, 3, 4};
//+
Physical Curve("CD", 6) = {3};
//+
Physical Point("A", 7) = {1};
//+
Physical Point("B", 8) = {2};
//+
Physical Point("C", 9) = {3};
//+
Physical Point("D", 10) = {4};
//+
Physical Point("E", 11) = {5};
