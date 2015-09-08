<h1> <span class="SC">buff-em</span> core library documentation
</h1>

The [[buff-em]] core library, `libbuff` exports a C++ class named 
`SWGGeometry.` The public methods of this class offer access to the 
computational routines provided by [[buff-em]] to implement the 
discretized volume-integral-equation (VIE) approach to 
computational electromagnetism.

The `SWGGeometry` class is analogous to the `RWGGeometry` class in 
the [<span class="SC">scuff-em</span> core library][libscuff]

# Creating an SWGGeometry

````C
  SWGGeometry *G = new SWGGeometry("E10Sphere_533.buffgeo");
````

# Setting up and solving the VIE system

Note: Incident fields are handled the same way as in [[scuff-em]]; see
[Incident fields in <span class="SC">scuff-em</span>][IncField].

````C
  HMatrix *M     = G->AllocateVIEMatrix();
  HVector *J     = G->AllocateRHSVector();

  cdouble E0[3]  = {1.0, 0.0, 0.0};
  double nHat[3] = {0.0, 0.0, 1.0};
  PlaneWave *PW  = new PlaneWave(E0, nHat);

  G->AssembleVIEMatrix(Omega, M);
  G->AssembleRHSVector(Omega, PW, J);
  M->LUFactorize();
  M->LUSolve(J);
````

# Computing scattered and total fields

````C

  HMatrix *XMatrix=new HMatrix("ListOfEvaluationPoints");

  HMatrix *EHIncMatrix = G->GetFields(PW,  J, Omega, XMatrix);
  HMatrix *EHScatMatrix = G->GetFields(0,  J, Omega, XMatrix);

  // E, H fields at the nth evaluation point are now in the 
  // nth row of the matrices, like this: Ex Ey Ez Hx Hy Hz
  
  int np=7;
  double X[3];
  cdouble E[3], H[3];

  XMatrix->GetEntries(np, ":", X);
  E[0] = EHIncMatrix->GetEntry(np, 0) + EHScatMatrix->GetEntry(np, 0);
  E[1] = EHIncMatrix->GetEntry(np, 1) + EHScatMatrix->GetEntry(np, 1);
  E[2] = EHIncMatrix->GetEntry(np, 2) + EHScatMatrix->GetEntry(np, 2);
  H[0] = EHIncMatrix->GetEntry(np, 3) + EHScatMatrix->GetEntry(np, 3);
  H[1] = EHIncMatrix->GetEntry(np, 4) + EHScatMatrix->GetEntry(np, 4);
  H[2] = EHIncMatrix->GetEntry(np, 5) + EHScatMatrix->GetEntry(np, 5);

  // now X[0..2] are the coordinates of evaluation point #np
  // and E[0..2], H[0..2] are the total E and H fields there
````

# Transforming the geometry 

To rotate the object named "SmallerTorus" by 90 degrees around the *z*
axis, then displaced it upward by 3 length units:
````C
 SWGVolume *O = G->GetObjectByLabel("SmallerTorus");
 if (!O) 
  ErrExit("no such object"); 

 O->Transform("ROTATED 90 ABOUT 0 0 1");
 O->Transform("DISPLACED 0 0 3");
````

# Computing power, force, and torque

````C
 PFTOptions MyOptions;
 InitPFTOptions(&MyOptions);

 PFTOptions->PFTMethod = BUFF_PFT_JDE; // or BUFF_PFT_DSI or BUFF_PFT_OVERLAP
 PFTOptions->IF = PW; // incident field is needed for some PFT methods

 HMatrix *PFTMatrix=G->GetPFT(J, Omega);

 // PFTMatrix->GetEntryD(no, nq) now returns the #nqth power quantity for object #no
 // no=0,1, ..., G->NumObjects - 1
 // nq=0,1   for absorbed, scattered power
 //    2,3,4 for force (x,y,z)
 //    5,6,7 for torque (x,y,z)

````

libscuff:	http://homerreid.github.io/scuff-em-documentation/API/libscuff/
IncField:	http://homerreid.com/scuff-EM/libscuff/IncField.shtml
