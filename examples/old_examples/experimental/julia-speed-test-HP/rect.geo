
Point(1) = {0, 20, 0, 20.0};
Point(2) = {0, 0, 0, 20.0};
Point(3) = {40, 0, 0, 20.0};
Point(4) = {40, 20, 0, 20.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};  

Plane Surface(1) = {1};

Transfinite Curve {1, 3} = 1001 Using Progression 1;
Transfinite Curve {2, 4} = 1001 Using Progression 1;

Transfinite Surface {1};

Recombine Surface {1};

Mesh.ElementOrder = 1;
Mesh.SecondOrderIncomplete = 1;

Mesh 2;
