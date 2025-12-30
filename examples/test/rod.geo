//+
SetFactory("OpenCASCADE");
//+
length = DefineNumber[ 100, Name "Parameters/length" ];
//+
elements = DefineNumber[ 1, Name "Parameters/elements" ];
//+
Rectangle(1) = {0, 0, 0, length, length/2, 0};
//+
//Transfinite Curve {1, 3} = elements+1 Using Progression 1;
//+
//Transfinite Curve {4, 2} = 2 Using Progression 1;
//+
//Transfinite Surface {1};
//+
//Recombine Surface {1};
//+
//Mesh.ElementOrder = 1;
//+
MeshSize {:} = elements;
Mesh 2;
//+
Physical Curve("bottom", 5) = {1};
//+
Physical Curve("left", 6) = {4};

//+
Physical Curve("right", 7) = {2};
//+
Physical Surface("body", 8) = {1};
//+
//Point(5) = {0, 0.5, 0, 1.0};
//+
//Point(6) = {length, 0.5, 0, 1.0};
//+
//Line(5) = {5, 6};
//+
//Physical Curve("path", 9) = {5};
