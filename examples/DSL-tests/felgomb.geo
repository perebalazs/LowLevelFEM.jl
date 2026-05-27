//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Point(1) = {-1, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Circle(1) = {1, 2, 3};
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; 
}
//+
MeshSize {1, 3, 4, 2} = 0.2;
Mesh.ElementOrder=2;
Mesh 2;
//+
Physical Surface("shell", 5) = {1};
//+
Physical Curve("x", 6) = {3};
//+
Physical Curve("y", 7) = {4};
//+
Physical Curve("z", 8) = {1};
