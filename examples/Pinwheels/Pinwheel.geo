//////////////////////////////////////////////////
// GMSH geometry for a QED pinwheel
// 
// Homer Reid 8/2015
//
// usage:
// 
// 1. To produce a volume mesh:
//  gmsh -3 [ options ] Pinwheel.geo
// 
// 2. To produce a surface mesh:
//  gmsh -2 -setnumber Mesh3D 0 [ options ] Pinwheel.geo
// 
// options:
//  --setnumber NumArms   2--7       (number of arms)
//  --setnumber Theta     Pi/6       (angular extent of arm)
//  --setnumber R1        0.4        (radius of inner disc)
//  --setnumber R2        0.2        (rounding radius for pinwheel arms)
//  --setnumber L1        1.0        (primary arm length)
//  --setnumber L2        0.2        (secondary arm length)
//  --setnumber T         0.2        (slab thickness)
//  --setnumber Mesh3D    1          (for volume mesh)
//////////////////////////////////////////////////

DefineConstant[ NumArms       = 3    ];
DefineConstant[ Theta         = Pi/6 ];
DefineConstant[ R1            = 0.4  ];
DefineConstant[ R2            = 0.2  ];
DefineConstant[ L1            = 1.0  ];
DefineConstant[ L2            = 0.2  ];
DefineConstant[ T             = 0.2  ];
DefineConstant[ Mesh3D        = 1    ];

l=0.2;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
X00=0.0;
Y00=0.0;

X10=R1*Cos(Theta);
Y10=R1*Sin(Theta);

X11=X10+L1;
Y11=Y10;

X01=X11;
Y01=Y11 + R2;

X12=X01 + R2;
Y12=Y01;

X13=X12;
Y13=Y12+L2;

X02=X13 + R2;
Y02=Y13;

X14=X02;
Y14=Y02 + R2;

X15=X02 + R2;
Y15=Y02;

X18=R1*Cos(Theta);
Y18=-R1*Sin(Theta);

X16=X15;
Y16=Y18 + R2;

X03=X16 - R2;
Y03=Y16;

X17=X03;
Y17=Y18;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Point(00)  = { X00, Y00,  0.5*T,  l};
Point(01)  = { X01, Y01,  0.5*T,  l};
Point(02)  = { X02, Y02,  0.5*T,  l};
Point(03)  = { X03, Y03,  0.5*T,  l};

Point(10)  = { X10, Y10,  0.5*T,  l};
Point(11)  = { X11, Y11,  0.5*T,  l};
Point(12)  = { X12, Y12,  0.5*T,  l};
Point(13)  = { X13, Y13,  0.5*T,  l};
Point(14)  = { X14, Y14,  0.5*T,  l};
Point(15)  = { X15, Y15,  0.5*T,  l};
Point(16)  = { X16, Y16,  0.5*T,  l};
Point(17)  = { X17, Y17,  0.5*T,  l};
Point(18)  = { X18, Y18,  0.5*T,  l};

Line(1)   = {10, 11};
Circle(2) = {11, 01, 12};
Line(3)   = {12, 13};
Circle(4) = {13, 02, 14};
Circle(5) = {14, 02, 15};
Line(6)   = {15, 16};
Circle(7) = {16, 03, 17};
Line(8)   = {17, 18};
Circle(9) = {18, 00, 10};

Line Loop(1)     = {1,2,3,4,5,6,7,8,9};
Ruled Surface(1) = {1};

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
Phi = 2*Pi/NumArms;
Point(19)          = {R1*Cos(Phi-Theta), R1*Sin(Phi-Theta),  0.5*T,  l};

Circle(10)         = {10, 00, 19};
Line(11)           = {00, 18};
Line(12)           = {19, 00};
Line Loop(2)       = {11, 9, 10, 12};
Ruled Surface(2)   = {2};

Extrude{0,0,-T} { Line{1,2,3,4,5,6,7,8,10}; }
Circle(49) = { 39, 41, 20 };

Line Loop(3) = {13, 17, 21, 25, 29, 33, 37, 41, 49};
Ruled Surface(3) = {3};

Circle(50)         = {20, 41, 42}; 
Line(51)           = {42, 41};
Line(52)           = {41, 39};
Line Loop(4)       = {49, 50, 51, 52};
Ruled Surface(4)   = {4};


//////////////////////////////////////////////////
//* instructions for producing a surface mesh 
//////////////////////////////////////////////////

If (Mesh3D == 0)

   For m In {1:(NumArms-1)}
    Rotate{ {0, 0, 1}, {0,0,0}, m*Phi } 
    { Duplicata{Surface{1,2,3,4,16,20,24,28,32,36,40,44,48};}}

  EndFor

EndIf

//////////////////////////////////////////////////
// instructions for producing a volume mesh: extrude
// the uppermost surfaces in the thickness direction
// to form volumes, then duplicate those by rotation
//////////////////////////////////////////////////
If (Mesh3D == 1)
  Extrude{0,0,-T} { Surface{1, 2}; }

  If( NumArms>1 )
    For m In {1:(NumArms-1)}
     Rotate{ {0, 0, 1}, {0,0,0}, m*Phi } { Duplicata{Volume{1,2}; } }
    EndFor
  EndIf

EndIf
