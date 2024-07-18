// Gmsh project created on Fri Sep 08 16:11:50 2023
SetFactory("OpenCASCADE");
//+

//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 2, 0, 1.0};
//+
Point(3) = {2, 2, 0, 1.0};
//+
Point(4) = {2, 0, 0, 1.0};
//+
Point(5) = {0, 0, 2, 1.0};
//+
Point(6) = {0, 2, 2, 1.0};
//+
Point(7) = {2, 2, 2, 1.0};
//+
Point(8) = {2, 0, 2, 1.0};
//+
Point(9) = {1, 1, 1, 1.0};
//+


//+
Point(10) = {1, 1, 0, 1.0};
//+
Point(11) = {2, 1, 1, 1.0};
//+
Point(12) = {1, 2, 1, 1.0};
//+
Point(13) = {0, 1, 1, 1.0};
//+
Point(14) = {1, 0, 1, 1.0};
//+
Point(15) = {1, 1, 2, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 5};
//+
Line(7) = {5, 6};
//+
Line(8) = {6, 7};
//+
Line(9) = {8, 4};
//+
Line(10) = {7, 3};
//+
Line(11) = {6, 2};
//+
Line(12) = {5, 1};
//+
Line(13) = {1, 10};
//+
Line(14) = {10, 4};
//+
Line(15) = {10, 2};
//+
Line(16) = {10, 3};
//+
Line(17) = {10, 9};
//+
Line(18) = {9, 1};
//+
Line(19) = {9, 4};
//+
Line(20) = {9, 3};
//+
Line(21) = {9, 2};
//+
Curve Loop(1) = {14, -1, 13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 15, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {16, 3, -15};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {14, 2, -16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, -19, -17};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {13, 17, 18};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {16, -20, -17};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {15, -21, -17};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {19, 2, -20};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {1, -19, 18};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {20, 3, -21};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {18, -4, -21};
//+
Plane Surface(12) = {12};
//+
Surface Loop(1) = {4, 9, 5, 7};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {1, 10, 6, 5};
//+
Surface Loop(3) = {1, 10, 6, 5};
//+
Volume(2) = {3};
//+
Surface Loop(4) = {2, 12, 8, 6};
//+
Volume(3) = {4};
//+
Surface Loop(5) = {3, 11, 7, 8};
//+
Volume(4) = {5};



Line(22) = {3, 12};

//+
Line(23) = {12, 6};
//+
Line(24) = {12, 2};
//+
Line(25) = {12, 7};
//+
Line(26) = {12, 9};
//+
Line(27) = {7, 9};
//+
Line(28) = {9, 6};

Transfinite Line {:}=2;


//+
Curve Loop(13) = {22, 25, 10};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {25, -8, -23};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {23, 11, -24};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {22, 24, -3};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {20, 22, 26};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {26, -27, -25};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {26, 28, -23};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {24, -21, -26};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {20, -10, 27};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {27, 28, 8};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {28, 11, -21};
//+
Plane Surface(23) = {23};
//+
Surface Loop(6) = {21, 13, 17, 18};
//+
Volume(5) = {6};
//+
Surface Loop(7) = {14, 22, 19, 18};
//+
Volume(6) = {7};
//+
Surface Loop(8) = {15, 23, 20, 19};
//+
Volume(7) = {8};
//+


//+
Surface Loop(10) = {16, 17, 20, 11};
//+
Volume(9) = {10};
//+
Line(29) = {13, 2};
//+
Line(30) = {13, 5};
//+
Line(31) = {13, 1};
//+
Line(32) = {13, 6};
//+
Line(33) = {13, 9};
//+
Curve Loop(24) = {29, -11, -32};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {32, -7, -30};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {31, -12, -30};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {29, 4, -31};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {21, -29, 33};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {32, -28, -33};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {33, 18, -31};
//+
Plane Surface(30) = {30};
//+
Line(34) = {5, 9};
//+
Curve Loop(31) = {33, -34, -30};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {28, -7, 34};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {34, 18, -12};
//+
Plane Surface(33) = {33};
//+

//+
Surface Loop(11) = {24, 23, 28, 29};
//+
Volume(10) = {11};



//+
Surface Loop(12) = {28, 27, 12, 30};
//+
Volume(11) = {12};
//+
Surface Loop(13) = {26, 33, 30, 31};
//+
Volume(12) = {13};
//+
Surface Loop(14) = {25, 32, 31, 29};
//+
Volume(13) = {14};
//+
Line(35) = {15, 7};
//+
Line(36) = {15, 8};
//+
Line(37) = {5, 15};
//+
Line(38) = {6, 15};
//+
Line(39) = {15, 9};
//+
Curve Loop(34) = {37, -38, -7};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {8, -35, -38};
//+
Plane Surface(35) = {35};
//+
Curve Loop(36) = {35, 5, -36};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {6, 37, 36};
//+
Plane Surface(37) = {37};
//+
Curve Loop(38) = {38, 39, 28};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {37, 39, -34};
//+
Plane Surface(39) = {39};
//+
Curve Loop(40) = {39, -27, -35};
//+
Plane Surface(40) = {40};
//+
Line(40) = {8, 9};
//+
Curve Loop(41) = {36, 40, -39};
//+
Plane Surface(41) = {41};
//+
Surface Loop(15) = {34, 32, 39, 38};
//+
Volume(14) = {15};
//+
Surface Loop(16) = {22, 35, 40, 38};
//+
Volume(15) = {16};
//+
Curve Loop(42) = {5, 40, -27};
//+
Plane Surface(42) = {42};
//+
Curve Loop(43) = {40, -34, -6};
//+
Plane Surface(43) = {43};
//+
Surface Loop(17) = {42, 36, 40, 41};
//+
Volume(16) = {17};
//+
Surface Loop(18) = {41, 37, 43, 39};
//+
Volume(17) = {18};




//+
Line(41) = {5, 14};
//+
Line(42) = {14, 4};
//+
Line(43) = {14, 8};
//+
Line(44) = {14, 1};
//+
Line(45) = {14, 9};
//+
Curve Loop(44) = {12, -44, -41};
//+
Plane Surface(44) = {44};
//+
Curve Loop(45) = {44, 1, -42};
//+
Plane Surface(45) = {45};
//+
Curve Loop(46) = {42, -9, -43};
//+
Plane Surface(46) = {46};
//+
Curve Loop(47) = {43, 6, 41};
//+
Plane Surface(47) = {47};
//+
Curve Loop(48) = {41, 45, -34};
//+
Plane Surface(48) = {48};
//+
Curve Loop(49) = {43, 40, -45};
//+
Plane Surface(49) = {49};
//+
Curve Loop(50) = {45, 19, -42};
//+
Plane Surface(50) = {50};
//+
Curve Loop(51) = {44, -18, -45};
//+
Plane Surface(51) = {51};
//+
Curve Loop(52) = {9, -19, -40};
//+
Plane Surface(52) = {52};
//+
Surface Loop(19) = {43, 47, 49, 48};
//+
Volume(18) = {19};
//+
Surface Loop(20) = {44, 33, 51, 48};
//+
Volume(19) = {20};
//+
Surface Loop(21) = {51, 45, 10, 50};
//+
Volume(20) = {21};



//+
Surface Loop(22) = {46, 52, 49, 50};
//+
Volume(21) = {22};
//+
Line(46) = {11, 7};
//+
Line(47) = {11, 3};
//+
Line(48) = {11, 4};
//+
Line(49) = {11, 8};
//+
Line(50) = {9, 11};
//+
Curve Loop(53) = {10, -47, 46};
//+
Plane Surface(53) = {53};
//+
Curve Loop(54) = {46, 5, -49};
//+
Curve Loop(55) = {5, -49, 46};
//+
Plane Surface(54) = {55};
//+
Curve Loop(56) = {9, -48, 49};
//+
Plane Surface(55) = {56};
//+
Curve Loop(57) = {48, 2, -47};
//+
Plane Surface(56) = {57};
//+
Curve Loop(58) = {46, 27, 50};
//+
Plane Surface(57) = {58};
//+
Curve Loop(59) = {50, 47, -20};
//+
Plane Surface(58) = {59};
//+
Curve Loop(60) = {19, -48, -50};
//+
Plane Surface(59) = {60};
//+
Curve Loop(61) = {50, 49, 40};
//+
Plane Surface(60) = {61};

Surface Loop(23) = {55, 52, 60, 59};
//+
Volume(22) = {23};
//+
Surface Loop(24) = {54, 42, 57, 60};
//+
Volume(23) = {24};
//+
Surface Loop(25) = {56, 9, 58, 59};
//+
Volume(24) = {25};
//+
Surface Loop(26) = {21, 53, 57, 58};
//+
Volume(25) = {26};



Transfinite Line {:}=2;

Physical Volume("body", 62) = {7, 4, 6, 15, 23, 16, 14, 13, 25, 5, 10, 17, 18, 9, 21, 22, 19, 12, 1, 24, 11, 3, 20, 2};

Mesh 3;