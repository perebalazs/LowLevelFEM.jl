//+ LowLevelFEM
//+
SetFactory("OpenCASCADE");
//+
b = DefineNumber[ 50, Name "Parameters/base" ];
//+
h = DefineNumber[ 100, Name "Parameters/height" ];
//+
t = DefineNumber[ 5, Name "Parameters/thikness" ];
//+
l = DefineNumber[ 500, Name "Parameters/length" ];
//+
n = DefineNumber[ 5, Name "Parameters/elements" ];
//+
p = DefineNumber[ 1, Name "Parameters/order" ];
//+
Point(1) = {b, h/2-t/2, 0, 1.0};
//+
Point(2) = {t/2, h/2-t/2, 0, 1.0};
//+
Point(3) = {t/2, -h/2+t/2, 0, 1.0};
//+
Point(4) = {b, -h/2+t/2, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Transfinite Curve {1, 3} = n+1 Using Progression 1;
n1 = Ceil((h-t) / (b / n));
Transfinite Curve {2} = n1+1 Using Progression 1;

Mesh 1;


//+
Extrude {0, 0, l} {
  Curve{1}; Curve{2}; Curve{3}; Layers {Ceil(l / ((b-t/2)/n)) + 1}; Recombine;
}

//Point(100) = {-11.5, 0, l, 1.0};

Mesh.ElementOrder=p;
Mesh 2;
//+
Physical Surface("shell", 17) = {1, 2, 3};
//+
Physical Curve("support", 18) = {1, 2, 3};
//+
Physical Curve("load", 19) = {6, 8, 10};
//+
//Physical Point("remote", 20) = {100};
