//+
Point(1) = {0, 100, 0, 1.0};
//+
Point(2) = {200, 100, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Line(1) = {3, 2};
//+
Line(2) = {1, 2};
//+
Transfinite Curve {2, 1} = 2 Using Progression 1;
//+
Mesh 1;
//+
Physical Curve("rod", 3) = {2, 1};
//+
Physical Point("A", 4) = {3};
//+
Physical Point("B", 5) = {1};
//+
Physical Point("C", 6) = {2};
