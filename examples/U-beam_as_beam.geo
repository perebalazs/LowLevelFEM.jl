//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
l = DefineNumber[ 500, Name "Parameters/length" ];
//+
n = DefineNumber[ 10, Name "Parameters/elements" ];
//+
p = DefineNumber[ 2, Name "Parameters/order" ];
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0, l, 1.0};
//+
Line(1) = {1, 2};

//+
Transfinite Curve {1} = n+1 Using Progression 1;

Mesh.ElementOrder=p;
Mesh 1;

//+
Physical Curve("beam", 2) = {1};
//+
Physical Point("support", 3) = {1};
//+
Physical Point("load", 4) = {2};
