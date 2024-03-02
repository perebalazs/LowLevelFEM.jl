SetFactory("OpenCASCADE");

a = 3.5;


//+
Box(1) = {0, 0, 0, a, a, a};

R = 0.705;

falvastagsag = 0.023;


// 1-es gömbhéj
Sphere(2) = {2.75,2.75,2.75, R,  -Pi/2, Pi/2, 2*Pi};
//+
Sphere(3) = {2.75,2.75,2.75,R-falvastagsag,  -Pi/2, Pi/2, 2*Pi};


// 2-es gömbhéj 
Sphere(4) = {0.72,2.75,2.75, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(5) = {0.72,2.75,2.75, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 3-mas gömbhéj
Sphere(6) = {2.75,0.72,2.75, R, -Pi/2, Pi/2, 2*Pi};  // R változtatva_!!!
//+
Sphere(7) = {2.75,0.72,2.75, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 4-es gömbhéj
Sphere(8) = {0.72,0.72,2.75, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(9) = {0.72,0.72,2.75, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};

// 5-ös gömbhéj
Sphere(10) = {1.73,0.75,1.74, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(11) = {1.73,0.75,1.74, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};

// 6-os gömbhéj
Sphere(12) = {1.75,2.75,1.74, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(13) = {1.75,2.75,1.74, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 7-es gömbhéj
Sphere(14) = {2.75,1.73,1.74, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(15) = {2.75,1.73,1.74, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 8-as gömbhéj
Sphere(16) = {0.72,1.75,1.74, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(17) = {0.72,1.75,1.74, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};





// 9-es gömbhéj
Sphere(18) = {2.75,2.75,0.73, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(19) = {2.75,2.75,0.73, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 10-es gömbhéj
Sphere(20) = {0.72,2.75,0.73, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(21) = {0.72,2.75,0.73, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 11-es gömbhéj
Sphere(22) = {2.75,0.72,0.73, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(23) = {2.75,0.72,0.73, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};

// 12-es gömbhéj
Sphere(24) = {0.72,0.72,0.73, R, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(25) = {0.72,0.72,0.73, R-falvastagsag, -Pi/2, Pi/2, 2*Pi};







// 13-mas gömbhéj
Sphere(26) = {1.75,1.78,2.75, R/1.06, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(27) = {1.75,1.78,2.75, R/1.06-falvastagsag, -Pi/2, Pi/2, 2*Pi};


// 14-es gömbhéj
Sphere(28) = {1.75,1.78,0.73, R/1.06, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(29) = {1.75,1.78,0.73, R/1.06-falvastagsag, -Pi/2, Pi/2, 2*Pi};



Delete { Volume{:}; }




// 1-től 6-ig megy itt a a kocka felület számozása
Surface Loop(30) = {6, 1, 3, 5, 2, 4};  // (1+gömbök száma+1 -től indul)


// az összes külső és belső gömb  felületét átalakítjuk 7-től megy felfelé az összes
// mindet át kell alakítani!
Surface Loop(31) = {7};
Surface Loop(32) = {8};

Surface Loop(33) = {9};
Surface Loop(34) = {10};

Surface Loop(35) = {11};
Surface Loop(36) = {12};

Surface Loop(37) = {13};
Surface Loop(38) = {14};

Surface Loop(39) = {15};
Surface Loop(40) = {16};

Surface Loop(41) = {17};
Surface Loop(42) = {18};

Surface Loop(43) = {19};
Surface Loop(44) = {20};

Surface Loop(45) = {21};
Surface Loop(46) = {22};

Surface Loop(47) = {23};
Surface Loop(48) = {24};

Surface Loop(49) = {25};
Surface Loop(50) = {26};

Surface Loop(51) = {27};
Surface Loop(52) = {28};

Surface Loop(53) = {29};
Surface Loop(54) = {30};

Surface Loop(55) = {31};
Surface Loop(56) = {32};

Surface Loop(57) = {33};
Surface Loop(58) = {34};




// gömbhélyak között páronként megcsináljuk a térfogatot

Volume(1) = {31,32};
Volume(2) = {33,34};
Volume(3) = {35,36};
Volume(4) = {37,38};
Volume(5) = {39,40};
Volume(6) = {41,42};
Volume(7) = {43,44};
Volume(8) = {45,46};
Volume(9) = {47,48};
Volume(10) = {49,50};
Volume(11) = {51,52};
Volume(12) = {53,54};
Volume(13) = {55,56};
Volume(14) = {57,58};



//Létre hozzuk a külső gömbhéj és a kocka közötti térfogatot. AZ UTOLSÓ A KOCKA
Volume(16) = {31,33,35,37,39,41,43,45,47,49,51,53,55,57,30};


Transfinite Surface {1,2,3,4,5,6};  // felület kisimítása a kockán
Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12}=10;  // kocka vonalai
Transfinite Line {14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95}=40;  // gömtöb vonalának felosztása  // 14-TŐL 3-MASÁVAL UGRÁL
Mesh 3;


Physical Volume(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};  // itt a {}-közé bepakolod a 2 gömb közti térfogatokat

Save "save_volume_1.m";