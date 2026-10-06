//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Box(2) = {0, 0, 1, 1, 1, 1};
//+
Coherence;
//+
Mesh 3;
SetOrder 5;
//+
Physical Volume("part1", 21) = {1};
//+
Physical Volume("part2", 22) = {2};
