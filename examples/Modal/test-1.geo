//+
SetFactory("OpenCASCADE");
Rectangle(1) = {100, 0, 0, 100, 1000, 0};
//+
Mesh.ElementOrder=1;
//+
Mesh 2;
//+
Physical Surface("body", 5) = {1};
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("right", 7) = {2};
//+
Physical Curve("top", 8) = {3};
//+
Physical Curve("bottom", 9) = {1};

