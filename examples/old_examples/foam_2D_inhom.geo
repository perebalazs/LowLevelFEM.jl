// Gmsh project created on Sat Aug 05 12:03:04 2023
SetFactory("OpenCASCADE");
//+

lc=1; 

Point(1) = {0, 0, 0, lc};
//+
Point(2) = {0, 2.5, 0, lc};
//+
Point(3) = {2.5, 2.5, 0, lc};
//+
Point(4) = {2.5, 0, 0, lc};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};


Curve Loop(1) = {1,2,3,4};  // vonalakból
//+
//Plane Surface(1) = {1};


//+
Circle(5) = {0.3, 0.29, 0, 0.1+0.15, 0, 2*Pi};
//+
Circle(6) = {0.3, 0.29, 0, 0.08+0.15, 0, 2*Pi};
//+

Circle(7) = {0.4, 0.85, 0, 0.11+0.15, 0, 2*Pi};
//+
Circle(8) = {0.4, 0.85, 0, 0.09+0.15, 0, 2*Pi};
//+

Circle(9) = {0.95, 0.87, 0, 0.12+0.15, 0, 2*Pi};
//+
Circle(10) = {0.95, 0.87, 0, 0.1+0.15, 0, 2*Pi};
//+

Circle(11) = {0.9, 0.3, 0, 0.12+0.15, 0, 2*Pi};
//+
Circle(12) = {0.9, 0.3, 0, 0.1+0.15, 0, 2*Pi};
//+

Circle(13) = {1.3, 1.35, 0, 0.12+0.15, 0, 2*Pi};
//+
Circle(14) = {1.3, 1.35, 0, 0.1+0.15, 0, 2*Pi};
//+

Circle(15) = {1.5, 0.85, 0, 0.1+0.15, 0, 2*Pi};
//+
Circle(16) = {1.5, 0.85, 0, 0.08+0.15, 0, 2*Pi};
//+

Circle(17) = {0.72, 1.36, 0, 0.1+0.15, 0, 2*Pi};
//+
Circle(18) = {0.72, 1.36, 0, 0.08+0.15, 0, 2*Pi};
//+

Circle(19) = {1.48, 0.32, 0, 0.12+0.12, 0, 2*Pi};
//+
Circle(20) = {1.48, 0.32, 0, 0.1+0.12, 0, 2*Pi};
//+


Circle(21) = {1.85, 1.3, 0, 0.11+0.15, 0, 2*Pi};
//+
Circle(22) = {1.85, 1.3, 0, 0.09+0.15, 0, 2*Pi};
//+


Circle(23) = {2.1, 0.7, 0, 0.12+0.2, 0, 2*Pi};
//+
Circle(24) = {2.1, 0.7, 0, 0.1+0.2, 0, 2*Pi};
//+

Circle(25) = {1.9, 0.2, 0, 0.12+0.05, 0, 2*Pi};
//+
Circle(26) = {1.9, 0.2, 0, 0.1+0.05, 0, 2*Pi};
//+

Circle(27) = {2.3, 0.2, 0, 0.12+0.05, 0, 2*Pi};
//+
Circle(28) = {2.3, 0.2, 0, 0.1+0.05, 0, 2*Pi};
//+


Circle(29) = {2.3, 1.4, 0, 0.12+0.05, 0, 2*Pi};
//+
Circle(30) = {2.3, 1.4, 0, 0.1+0.05, 0, 2*Pi};
//+

Circle(31) = {0.26, 1.3, 0, 0.12+0.05, 0, 2*Pi};
//+
Circle(32) = {0.26, 1.3, 0, 0.1+0.05, 0, 2*Pi};
//+

Circle(33) = {0.3, 1.77, 0, 0.12+0.15, 0, 2*Pi};
//+
Circle(34) = {0.3, 1.77, 0, 0.1+0.15, 0, 2*Pi};
//+

Circle(35) = {2.2, 1.87, 0, 0.12+0.15, 0, 2*Pi};
//+
Circle(36) = {2.2, 1.87, 0, 0.1+0.15, 0, 2*Pi};
//+

Circle(37) = {0.87, 1.93, 0, 0.14+0.15, 0, 2*Pi};
//+
Circle(38) = {0.87, 1.93, 0, 0.12+0.15, 0, 2*Pi};
//+


Circle(39) = {1.55, 1.95, 0, 0.14+0.2, 0, 2*Pi};
//+
Circle(40) = {1.55, 1.95, 0, 0.12+0.2, 0, 2*Pi};
//+

Circle(41) = {0.3, 2.27, 0, 0.11+0.10, 0, 2*Pi};
//+
Circle(42) = {0.3, 2.27, 0, 0.09+0.10, 0, 2*Pi};
//+

Circle(43) = {1.15, 2.31, 0, 0.11+0.05, 0, 2*Pi};
//+
Circle(44) = {1.15, 2.31, 0, 0.09+0.05, 0, 2*Pi};

Circle(45) = {1.98, 2.29, 0, 0.11+0.07, 0, 2*Pi};
//+
Circle(46) = {1.98, 2.29, 0, 0.09+0.07, 0, 2*Pi};



Curve Loop(2) = {5};
Curve Loop(3) = {6};
Curve Loop(4) = {7};
Curve Loop(5) = {8};
Curve Loop(6) = {9};
Curve Loop(7) = {10};
Curve Loop(8) = {11};
Curve Loop(9) = {12};
Curve Loop(10) = {13};
Curve Loop(11) = {14};
Curve Loop(12) = {15};
Curve Loop(13) = {16};
Curve Loop(14) = {17};
Curve Loop(15) = {18};
Curve Loop(16) = {19};
Curve Loop(17) = {20};
Curve Loop(18) = {21};
Curve Loop(19) = {22};
Curve Loop(20) = {23};
Curve Loop(21) = {24};
Curve Loop(22) = {25};
Curve Loop(23) = {26};
Curve Loop(24) = {27};
Curve Loop(25) = {28};
Curve Loop(26) = {29};
Curve Loop(27) = {30};
Curve Loop(28) = {31};
Curve Loop(29) = {32};
Curve Loop(30) = {33};
Curve Loop(31) = {34};
Curve Loop(32) = {35};
Curve Loop(33) = {36};
Curve Loop(34) = {37};
Curve Loop(35) = {38};
Curve Loop(36) = {39};
Curve Loop(37) = {40};
Curve Loop(38) = {41};
Curve Loop(39) = {42};
Curve Loop(40) = {43};
Curve Loop(41) = {44};
Curve Loop(42) = {45};
Curve Loop(43) = {46};




Plane Surface(1) = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,1}; //nagy_buborek-negyzet

Plane Surface(2) = {2,3};			//buborek-buborek
Plane Surface(3) = {4,5};
Plane Surface(4) = {6,7};
Plane Surface(5) = {8,9};
Plane Surface(6) = {10,11};
Plane Surface(7) = {12,13};
Plane Surface(8) = {14,15};
Plane Surface(9) = {16,17};
Plane Surface(10) = {18,19};
Plane Surface(11) = {20,21};
Plane Surface(12) = {22,23};
Plane Surface(13) = {24,25};
Plane Surface(14) = {26,27};
Plane Surface(15) = {28,29};
Plane Surface(16) = {30,31};
Plane Surface(17) = {32,33};
Plane Surface(18) = {34,35};
Plane Surface(19) = {36,37};
Plane Surface(20) = {38,39};
Plane Surface(21) = {40,41};
Plane Surface(22) = {42,43};




Transfinite Curve {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46} = 100 Using Progression 1;
Transfinite Curve {1, 4, 3, 2} = 101 Using Progression 1;
Recombine Surface {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};

Mesh 2;

Physical Surface("Fe", 1) = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};  // itt a {}-közé bepakolod a 2 gömb közti térfogatokat




// 21-21-el náha jó
//+
Physical Surface("Al", 47) = {1};
//+
Physical Curve("supp", 48) = {4};
//+
Physical Curve("load", 49) = {2};
