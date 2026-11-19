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
Point(4) = {300, 100, 0, 1.0};
//+
Point(5) = {300, 0, 0, 1.0};
//+
Circle(1) = {1, 2, 3};
//+
Line(2) = {3, 4};
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; Curve{2}; 
}

//+
MeshSize {:} = 10;
Mesh.ElementOrder=2;
Mesh 2;
//+
Physical Surface("shell", 8) = {1, 2};
//+
Physical Curve("xy", 9) = {2, 1};
//+
Physical Curve("yz", 10) = {6};
//+
Physical Curve("zx", 11) = {7, 5};
