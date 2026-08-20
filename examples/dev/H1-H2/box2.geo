//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Wedge(1) = {-0.5, 0.2, -0.2, 0.5, 0.5, 0.5, 0};
//+
Physical Volume("body", 10) = {1};
//+
Physical Surface("left", 11) = {1};
//+
Physical Surface("right", 12) = {2};
