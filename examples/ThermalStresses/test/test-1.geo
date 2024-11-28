//+
SetFactory("OpenCASCADE");
Rectangle(1) = {100, 0, 0, 100, 100, 0};
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