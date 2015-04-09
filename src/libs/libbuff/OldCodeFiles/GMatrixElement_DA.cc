/***************************************************************/
/* get the dipole and quadupole moments of the current         */
/* distribution described by a single SWG basis function.      */
/***************************************************************/
void GetDQMoments(SWGVolume *O, int nf, double J[3], double Q[3][3],
                  bool NeedQ)
{
  SWGFace *F = O->Faces[nf];
  double A= F->Area;

  double *QP = O->Vertices + 3*F->iQP;
  double *QM = O->Vertices + 3*F->iQM;

  J[0] = 0.25*A*(QM[0] - QP[0]);
  J[1] = 0.25*A*(QM[1] - QP[1]);
  J[2] = 0.25*A*(QM[2] - QP[2]);

  if (!NeedQ) return;

  double PreFac = A/20.0;
  double *x0 = F->Centroid;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    Q[Mu][Nu] = PreFac * (   QM[Mu]*(QM[Nu]-x0[Nu])
                            -QP[Mu]*(QP[Nu]-x0[Nu])
                            +x0[Mu]*(QP[Nu]-QM[Nu])
                         );

}
/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions in the dipole approximation retaining NumTerms    */
/* terms (here NumTerms may be 1 or 2).                        */
/***************************************************************/
cdouble GetGMatrixElement_DA(SWGVolume *OA, int nfA, 
                             SWGVolume *OB, int nfB,
                             cdouble Omega, int NumTerms)
{

  SWGFace *FA = OA->Faces[nfA];
  SWGFace *FB = OB->Faces[nfB];

  /***************************************************************/
  /* get dipole moments ******************************************/
  /***************************************************************/
  double JA[3], JB[3];
  double QA[3][3], QB[3][3];
  cdouble GMuNu[3][3], GMuNuRho[3][3][3];

  if (NumTerms==1)
   { 
     GetDQMoments(OA, nfA, JA, QA, false);
     GetDQMoments(OB, nfB, JB, QB, false);
     CalcGC(FA->Centroid, FB->Centroid, Omega, 1.0, 1.0, GMuNu, 0, 0, 0);
   }
  else
   { GetDQMoments(OA, nfA, JA, QA, true);
     GetDQMoments(OB, nfB, JB, QB, true);
     CalcGC(FA->Centroid, FB->Centroid, Omega, 1.0, 1.0, GMuNu, 0, GMuNuRho, 0);
   };
  
  cdouble G=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    G += JA[Mu]*GMuNu[Mu][Nu]*JB[Nu];

  if (NumTerms==1) return G;

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    for(int Rho=0; Rho<3; Rho++)
     G += GMuNuRho[Mu][Nu][Rho]
           *( QA[Mu][Rho]*JB[Nu] - JA[Mu]*QB[Nu][Rho] );

  return G;
}
