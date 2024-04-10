SetFactory("OpenCASCADE");




//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {15, 0, 0, 1.0};
//+
Point(4) = {15, 5, 0, 1.0};
//+
Point(5) = {10, 5, 0, 1.0};
//+
Point(6) = {0, 5, 0, 1.0};

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
Line(6) = {6, 1};
//+
Line(7) = {2, 5};

//+
Circle(8) = {12.5, 2.5, 0, 1.2, 0, 2*Pi};
//+
Circle(9) = {12.5, 2.5, 0, 0.8, 0, 2*Pi};

//+
Curve Loop(1) = {8};
//+
Curve Loop(2) = {9};

//+
Curve Loop(3) = {1,7,5,6};  // vonalakból
Curve Loop(4) = {2,3,4,7};  // vonalakból


Plane Surface(2) = {1, 4};

Plane Surface(3) = {1,2};

Plane Surface(4) = {3};




//+
Transfinite Curve {2,3,4,7} = 12 Using Progression 1;
Transfinite Curve {1,5} = 10 Using Progression 1;
Transfinite Curve {6} = 10/2 Using Progression 1; // itt midig fele annyi elem legyen, mint feljebb!!!

//+
Transfinite Curve {8} = 26 Using Progression 1; // körök felosztása

//+
Transfinite Curve {9} = 26 Using Progression 1;
Recombine Surface {3};

Mesh.ElementOrder = 1;

Mesh 2;



// EGY ANYAG ESETEN:
//Physical Surface("body", 5) = {2,3,4};
//+
//Physical Curve("left", 6) = {6};
//+
//Physical Curve("right", 7) = {3};
//+
//Physical Curve("mid", 8) = {7};





// KET ANYAG ESETEN:
Physical Surface("Fe", 5) = {3};
//+
Physical Surface("Al", 6) = {2,4};
//+

Physical Curve("left", 7) = {6};
//+
Physical Curve("right", 8) = {3};
//+
Physical Curve("mid", 9) = {7};