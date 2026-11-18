//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 100, 0, 1.0};
//+
Point(4) = {500, 100, 0, 1.0};
//+
Point(5) = {500, 0, 0, 1.0};
//+
Point(6) = {600, 0, 0, 1.0};
//+
Circle(1) = {1, 2, 3};
//+
Circle(2) = {4, 5, 6};
//+
Line(3) = {3, 4};
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; Curve{3}; Curve{2}; 
}

//+
MeshSize {:} = 20;
Mesh 2;
//+
Physical Surface("shell", 11) = {1, 2, 3};
//+
Physical Curve("xy", 12) = {1, 3, 2};
//+
Physical Curve("zx", 13) = {6, 8, 10};
