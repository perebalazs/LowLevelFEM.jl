//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
x1 = 0; y1 = 0;
x2 = 1; y2 = 1;
x3 = 0; y3 = 1;
eps = Sqrt((x1-x2)^2+(y1-y2)^2) / 10000;
//+
Point(1) = {0, y1, 0, 1.0};
//+
Point(2) = {x2, y2, 0, 1.0};
//+
Point(3) = {x3, y3, 0, 1.0};
//+
Point(4) = {x3+eps, y3+eps, 0, 1.0};
//+
Line(1) = {1, 2};
Line(2) = {3, 4};
//+
Transfinite Curve {1} = 2 Using Progression 1;
Transfinite Curve {2} = 2 Using Progression 1;
//+
Mesh 1;
//+
Physical Curve("segment", 2) = {1};
//+
Physical Curve("node", 3) = {2};
//+
Physical Point("left", 4) = {1};
//+
Physical Point("right", 5) = {2};
