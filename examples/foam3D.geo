SetFactory("OpenCASCADE");




a = 2;


//+
Box(1) = {0, 0, 0, a, a, a};


// 1-es gömbhéj
Sphere(2) = {1,0.5,1, 0.44, -Pi/2, Pi/2, 2*Pi};
Dilate {{1,0.5,1}, {1, 1, 1.5}} { Volume{2}; }
Rotate {{1,0,0}, {1,0.5,1}, -Pi/7} { Volume{2}; }
//+
Sphere(3) = {1,0.5,1, 0.4, -Pi/2, Pi/2, 2*Pi};
Dilate {{1,0.5,1}, {1, 1, 1.5}} { Volume{3}; }
Rotate {{1,0,0}, {1,0.5,1}, -Pi/7} { Volume{3}; }



// 2-es gömbhéj 
Sphere(4) = {1,1.5,1, 0.44, -Pi/2, Pi/2, 2*Pi};
Dilate {{1,1.5,1}, {1.5, 1, 0.5}} { Volume{4}; }
Rotate {{0,1,0}, {1,1.5,1}, -Pi/7} { Volume{4}; }
//+
Sphere(5) = {1,1.5,1, 0.4, -Pi/2, Pi/2, 2*Pi};
Dilate {{1,1.5,1}, {1.5, 1, 0.5}} { Volume{5}; }
Rotate {{0,1,0}, {1,1.5,1}, -Pi/7} { Volume{5}; }


Delete { Volume{:}; }

// 1-től 6-ig megy itt a a kocka felület számozása
Surface Loop(6) = {6, 1, 3, 5, 2, 4};  // (1+gömbök száma+1 -től indul)


// az összes külső és belső gömb  felületét átalakítjuk 7-től megy felfelé az összes
// mindet át kell alakítani!
Surface Loop(7) = {7};
Surface Loop(8) = {8};

Surface Loop(9) = {9};
Surface Loop(10) = {10};


// gömbhélyak között páronként megcsináljuk a térfogatot

Volume(1) = {7,8};
Volume(2) = {9,10};

//Létre hozzuk a külső gömbhéj és a kocka közötti térfogatot. AZ UTOLSÓ A KOCKA
Volume(3) = {7,9,6};

Transfinite Surface {1,2,3,4,5,6};  // felület kisimítása a kockán
Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12}=10;  // kocka vonalai

Transfinite Line {14,17,20,23}=20;  // gömtöb vonalának felosztása  // 14-TŐL 3-MASÁVAL UGRÁL

Physical Volume(1) = {1,2};  // itt a {}-közé bepakolod a 2 gömb közti térfogatokat

//Save "save_volume_1.m";//+

SetName "foam";
//Mesh.ElementOrder = 1;
Mesh 3;
//Mesh.SaveAll=1;
//RenumberMeshNodes;
Physical Surface("supp", 25) = {1};
//+
Physical Surface("load", 26) = {2};
//+
Physical Volume("body", 27) = {3, 1, 2};
