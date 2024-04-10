//+
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {100, 0, 0, 1.0};
//+
Point(3) = {100, 49, 0, 1.0};
//+
Point(4) = {100, 51, 0, 1.0};
//+
Point(5) = {100, 100, 0, 1.0};
//+
Point(6) = {0, 100, 0, 1.0};
//+
Point(7) = {0, 50, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};

//+
Transfinite Curve {6, 7} = 51 Using Progression 1;
//+
Transfinite Curve {2, 4} = 50 Using Progression 1;
//+
Transfinite Curve {3} = 3 Using Progression 1;
//+
Transfinite Curve {1, 5} = 101 Using Progression 1;
//+
Transfinite Surface {1} = {1,2,5,6};
//+
Recombine Surface {1};


Mesh 2;


//+
Physical Surface("body", 7) = {1};
//+
Physical Curve("left", 8) = {6, 7};
//+
Physical Curve("source", 9) = {3};
//+
Physical Point("middle", 10) = {7};
//+
Point(8) = {100, 50, 0, 1.0};
//+
Line(8) = {7, 8};
//+
Physical Curve("path", 11) = {8};
