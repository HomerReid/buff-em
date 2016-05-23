//************************************************************
//* input parameters      
//************************************************************
DefineConstant[ R      = 1.0  ]; // radius of dielectric sphere
DefineConstant[ T      = 0.1  ]; // thickness of metal coating
DefineConstant[ Theta  = 90.0 ]; // deficit angle in degrees
DefineConstant[ Mesh3D = 1 ];    // =1 for volume mesh

//************************************************************
//* derived parameters 
//************************************************************
TT=Theta*Pi/180.0; // theta in radians
RR=R+T;            // outer radius

//************************************************************
//* meshing finenesses ***************************************
//************************************************************
lc = R/4;  // fineness at north pole

//************************************************************
//************************************************************
//************************************************************
Point(0) = { 0, 0, 0,   l};
Point(1) = { 0, 0, +R,  l};
Point(2) = {  R*Sin[TT], 0,  R*Cos[TT], 0.5*l};
Point(3) = { RR*Sin[TT], 0, RR*Cos[TT], 0.5*l};
Point(4) = { 0, 0, -RR, 0.5*l};

Circle(1) = { 1, 0, 2 };
Line  (2) = { 2, 3    };
Circle(3) = { 3, 0, 4 };
Line  (4) = { 4, 0, 1 };


If (Mesh3D==0)

  Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{1,2,3,4}; }
  Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{5,8,12,15}; }
  Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Line{18,21,25,28}; }

  // metal-air interface 
  Physical Surface(1) = {11, 24, 37, 7, 20, 33};

  // dielectric-air interface 
  Physical Surface(2) = {14, 27, 40};

  // dielectric-metal interface 
  Physical Surface(3) = {17, 30, 43};
EndIf

If (Mesh3D==1)

  Line Loop(1) = {1, 2, 3, 4};
  Ruled Surface(1) = {1};
  Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{1}; }

  Rotate{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Duplicata{Volume{1};} }
  Rotate{ {0,0,1}, {0,0,0}, 4*Pi/3 } { Duplicata{Volume{1};} }
  
EndIf
