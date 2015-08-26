/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of BUFF-EM.
 *
 * BUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * BUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * GetNEQMoments.cc 
 *
 * homer reid  -- 8/2015
 *
 */

#include <libscuff.h>
#include "buff-neq.h"

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Get1BFMoments(SWGVolume *O, int nf,
                   double JMu[3], double JMuNu[3][3])
{
  SWGFace *F  = O->Faces[nf];
  double  A   = F->Area;
  double *QP  = O->Vertices + 3*(F->iQP);
  double *QM  = O->Vertices + 3*(F->iQM);
  
  for(int Mu=0; Mu<3; Mu++)
   JMu[Mu] = A*(QM[Mu] - QP[Mu]) / 4.0;

// note: a more compact formula for the following would be nice
  double *V1  = O->Vertices + 3*(F->iV1);
  double *V2  = O->Vertices + 3*(F->iV2);
  double *V3  = O->Vertices + 3*(F->iV3);
  double L1P[3], L2P[3], L3P[3], L1M[3], L2M[3], L3M[3];
  VecSub(V1, QP, L1P);
  VecSub(V2, QP, L2P);
  VecSub(V3, QP, L3P);
  VecSub(V1, QM, L1M);
  VecSub(V2, QM, L2M);
  VecSub(V3, QM, L3M);

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    JMuNu[Mu][Nu] 
     = + A*(  L1P[Mu]*L1P[Nu] - L1M[Mu]*L1M[Nu]
            + L2P[Mu]*L2P[Nu] - L2M[Mu]*L2M[Nu]
            + L3P[Mu]*L3P[Nu] - L3M[Mu]*L3M[Nu]
           ) / 30.0
       + A*(  L1P[Mu]*L2P[Nu] - L1M[Mu]*L2M[Nu]
             +L1P[Mu]*L3P[Nu] - L1M[Mu]*L3M[Nu]
             +L2P[Mu]*L3P[Nu] - L2M[Mu]*L3M[Nu]
             +L1P[Nu]*L2P[Mu] - L1M[Nu]*L2M[Mu]
             +L1P[Nu]*L3P[Mu] - L1M[Nu]*L3M[Mu]
             +L2P[Nu]*L3P[Mu] - L2M[Nu]*L3M[Mu]
           ) /60.0
       + A*(  L1P[Mu]*QP[Nu]  - L1M[Mu]*QM[Nu]
            + L2P[Mu]*QP[Nu]  - L2M[Mu]*QM[Nu]
            + L3P[Mu]*QP[Nu]  - L3M[Mu]*QM[Nu]
           ) /12.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetNEQMoments(BNEQData *BNEQD, int no, HMatrix *Rytov,
                   cdouble MMuNu[3][3], cdouble MMuNuRho[3][3][3])
{
  SWGGeometry *G    = BNEQD->G;
  SWGVolume *O      = G->Objects[no];
  //int Offset        = G->BFIndexOffset[no];
  int N             = O->NumInteriorFaces;
  HMatrix *SInverse = BNEQD->SInverse;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HVector *vMu[3], *vMuNu[3][3], *vTemp;
  cdouble *DataBuffer = BNEQD->WorkMatrix[1]->ZM;
  int nVec=0;
  for(int Mu=0; Mu<3; Mu++, nVec++)
   vMu[Mu] = new HVector(N, LHM_COMPLEX, DataBuffer + nVec*N);
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++, nVec++)
    vMuNu[Mu][Nu] = new HVector(N, LHM_COMPLEX, DataBuffer + nVec*N);
  vTemp=new HVector(N, LHM_COMPLEX, DataBuffer+nVec*N); nVec++;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("GNM Getting 1BF moments...");
  for(int n=0; n<N; n++)
   { 
     double JMu[3], JMuNu[3][3];
     Get1BFMoments(O, n, JMu, JMuNu);

     for(int Mu=0; Mu<3; Mu++)
      vMu[Mu]->SetEntry(n, JMu[Mu]);

     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       vMuNu[Mu][Nu]->SetEntry(n, JMuNu[Mu][Nu]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("GNM doing linear algebra...");
/*
  for(int Mu=0; Mu<3; Mu++)
   { vTemp->Copy(vMu[Mu]);
     SInverse->Apply(vTemp, vMu[Mu]);
   };

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { vTemp->Copy(vMuNu[Mu][Nu]);
      SInverse->Apply(vTemp, vMuNu[Mu][Nu]);
    };
*/

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    MMuNu[Mu][Nu]=Rytov->BilinearProduct(vMu[Mu], vMu[Nu]);

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    for(int Rho=0; Rho<3; Rho++)
     MMuNuRho[Mu][Nu][Rho]=Rytov->BilinearProduct(vMuNu[Mu][Rho], vMu[Nu]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int Mu=0; Mu<3; Mu++)
   delete vMu[Mu];
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    delete vMuNu[Mu][Nu];
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetMomentPFT(BNEQData *BNEQD, int no, double Omega,
                  HMatrix *Rytov, FILE *f)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble MMuNu[3][3], MMuNuRho[3][3][3];
  Log("GMP getting NEQ moments...");
  GetNEQMoments(BNEQD, no, Rytov, MMuNu, MMuNuRho);
  Log("GMP done...");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double F1[3], F2[3], F3[3], Torque[3];
  for(int Mu=0; Mu<3; Mu++)
   F1[Mu]=F2[Mu]=F3[Mu]=Torque[Mu]=0.0;

  double FPF = TENTHIRDS*ZVAC*Omega*Omega*Omega/(120.0*M_PI);
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { 
      F1[Mu] +=     FPF*imag( MMuNuRho[Nu][Mu][Nu] );
      F2[Mu] -= 3.0*FPF*imag( MMuNuRho[Mu][Nu][Nu] + MMuNuRho[Nu][Nu][Mu]);
      F3[Mu] += 5.0*FPF*imag( MMuNuRho[Mu][Nu][Nu] - MMuNuRho[Nu][Nu][Mu]);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double TPF = TENTHIRDS*ZVAC*Omega/(12.0*M_PI);
  for(int Mu=0; Mu<3; Mu++)
   { 
     int Nu=(Mu+1)%3, Rho=(Mu+2)%3;
     Torque[Mu] += TPF*imag( MMuNu[Nu][Rho] - MMuNu[Rho][Nu] );
   };


  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#if 0
  FILE *f=vfopen("%s.jl","a",GetFileBase(BNEQD->G->GeoFileName));
  static int Count=1;
  fprintf(f,"Omega_%i=%g;\n",Count,Omega);
  fprintf(f,"JJ_%i=im*zeros(3,3);\n",Count);
  fprintf(f,"MM_%i=im*zeros(3,3);\n",Count);
  fprintf(f,"Force_%i=%e;\n",Count,(F1[2]+F2[2]+F3[2])/FPF);
  fprintf(f,"Torque_%i=%e;\n",Count,Torque[2]/TPF);
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { 
      cdouble pp=Omega*Omega*MMuNu[Mu][Nu];
      fprintf(f,"pp_%i[%i,%i]=%e + (%e)*im;\n",
                 Count,Mu,Nu,real(pp),imag(pp));

      int Rho=(Nu+1)%3, Sigma=(Nu+2)%3;
      cdouble mp=II*Omega*(   MMuNuRho[Mu][Rho][Sigma]
                            - MMuNuRho[Mu][Sigma][Rho]
                          );
      fprintf(f,"mp_%i[%i,%i]=%e + (%e)*im;\n",
                 Count,Mu,Nu,real(mp),imag(mp));
    };
  fprintf(f,"\n");
  fclose(f);
  Count++;
#endif

  static HMatrix *ppMatrix = new HMatrix(3,3,LHM_COMPLEX);
  static HMatrix *mpMatrix = new HMatrix(3,3,LHM_COMPLEX);
  static HMatrix *U        = new HMatrix(3,3,LHM_COMPLEX);
  static HMatrix *VT       = new HMatrix(3,3,LHM_COMPLEX);
  static HVector *Lambda   = new HVector(3,LHM_REAL);
  static HVector *Sigma    = new HVector(3,LHM_REAL);

/*
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    ppMatrix->SetEntry(Mu,Nu, Omega*Omega*MMuNu[Mu][Nu]);
*/
  for(int Mu=0; Mu<3; Mu++)
   { ppMatrix->SetEntry(Mu, Mu, Omega*Omega*real(MMuNu[Mu][Mu]));
     for(int Nu=Mu+1; Nu<3; Nu++)
      { ppMatrix->SetEntry(Mu, Nu, Omega*Omega*MMuNu[Mu][Nu]);
        ppMatrix->SetEntry(Nu, Mu, Omega*Omega*conj(MMuNu[Mu][Nu]));
      };
   };

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { int NP1=(Nu+1)%3, NP2=(Nu+2)%3;
      mpMatrix->SetEntry(Mu,Nu, II*Omega*(  MMuNuRho[Mu][NP1][NP2]
                                          - MMuNuRho[Mu][NP2][NP1]
                                         )
                        );
    };

  cdouble p[3][3];
  ppMatrix->Eig(Lambda, U);
  double MaxEig=fmax( abs(Lambda->GetEntry(0)), abs(Lambda->GetEntry(1)));
  MaxEig=fmax( MaxEig, abs(Lambda->GetEntry(2)));
  for(int Mu=0.0; Mu<3.0; Mu++)
   if( abs(Lambda->GetEntry(Mu)) < 1.0e-8*MaxEig )
    Lambda->SetEntry(Mu, 0.0);

  for(int a=0; a<3; a++)
   for(int Mu=0; Mu<3; Mu++)
    p[a][Mu] = sqrt(Lambda->GetEntry(a))*U->GetEntry(Mu,a);

  cdouble m[3][3];
  mpMatrix->SVD(Sigma, U, VT);
  for(int a=0; a<3; a++)
   { 
     cdouble DotProd = conj(VT->GetEntry(a,0)) * p[a][0]
                      +conj(VT->GetEntry(a,1)) * p[a][1]
                      +conj(VT->GetEntry(a,2)) * p[a][2];
     cdouble LambdaA=Lambda->GetEntry(a);
     cdouble ScaleFactor 
      = LambdaA==0.0 ? 0.0 : Sigma->GetEntry(a) * DotProd / LambdaA;
     m[a][0] = ScaleFactor*DotProd*U->GetEntry(0,a);
     m[a][1] = ScaleFactor*DotProd*U->GetEntry(1,a);
     m[a][2] = ScaleFactor*DotProd*U->GetEntry(2,a);
   };

  double Fmxp[3]={0.0, 0.0, 0.0};
  double Tpxp[3]={0.0, 0.0, 0.0};
  TPF = TENTHIRDS*ZVAC/(12.0*M_PI*Omega);
  FPF = TPF*Omega;
  for(int a=0; a<3; a++)
   for(int Mu=0; Mu<3; Mu++)
    { int Nu=(Mu+1)%3, Rho=(Mu+2)%3;
      Fmxp[Mu]+=FPF*real( conj(m[a][Nu])*p[a][Rho] - conj(m[a][Rho])*p[a][Nu] );
      Tpxp[Mu]+=TPF*imag( conj(p[a][Nu])*p[a][Rho] - conj(p[a][Rho])*p[a][Nu] );
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
/*
  FILE *f=vfopen("%s.MomentPFT","r",BNEQD->FileBase);
  if (!f)
   { f=vfopen("%s.MomentPFT","w",BNEQD->FileBase);
     fprintf(f,"# 01       omega \n");
     fprintf(f,"# 02 03 04 Fx,Fy,Fz (term 1)\n");
     fprintf(f,"# 05 06 07 Fx,Fy,Fz (term 2)\n");
     fprintf(f,"# 08 09 10 Fx,Fy,Fz (term 3)\n");
     fprintf(f,"# 11 12 13 Tx,Ty,Tz (term 3)\n");
     fclose(f);
   };
  f=vfopen("%s.MomentPFT","a",BNEQD->FileBase);
  fprintf(f,"%e ",Omega);
  fprintf(f,"%e %e %e ",F1[0],F1[1],F1[2]);
  fprintf(f,"%e %e %e ",F2[0],F2[1],F2[2]);
  fprintf(f,"%e %e %e ",F3[0],F3[1],F3[2]);
  fprintf(f,"%e %e %e ",Torque[0],Torque[1],Torque[2]);
  fprintf(f,"%e %e %e ",Fmxp[0], Fmxp[1], Fmxp[2]);
  fprintf(f,"%e %e %e ",Tpxp[0], Tpxp[1], Tpxp[2]);
  fprintf(f,"\n");
  fclose(f);

  f=vfopen("%s.NEQMoments","a",BNEQD->FileBase);
  fprintf(f,"%e ",Omega);
  fprintf(f,"%e %e %e ",real(p[0][0]),real(p[0][1]),real(p[0][2]));
  fprintf(f,"%e %e %e ",imag(p[0][0]),imag(p[0][1]),imag(p[0][2]));
  fprintf(f,"%e %e %e ",real(p[1][0]),real(p[1][1]),real(p[1][2]));
  fprintf(f,"%e %e %e ",imag(p[1][0]),imag(p[1][1]),imag(p[1][2]));
  fprintf(f,"%e %e %e ",real(p[2][0]),real(p[2][1]),real(p[2][2]));
  fprintf(f,"%e %e %e ",imag(p[2][0]),imag(p[2][1]),imag(p[2][2]));
  fprintf(f,"%e %e %e ",real(m[0][0]),real(m[0][1]),real(m[0][2]));
  fprintf(f,"%e %e %e ",imag(m[0][0]),imag(m[0][1]),imag(m[0][2]));
  fprintf(f,"%e %e %e ",real(m[1][0]),real(m[1][1]),real(m[1][2]));
  fprintf(f,"%e %e %e ",imag(m[1][0]),imag(m[1][1]),imag(m[1][2]));
  fprintf(f,"%e %e %e ",real(m[2][0]),real(m[2][1]),real(m[2][2]));
  fprintf(f,"%e %e %e ",imag(m[2][0]),imag(m[2][1]),imag(m[2][2]));
  fclose(f);
*/
double PAbs1=0.0, PAbs2=0.0;
double PPF = ZVAC*Omega*Omega/(12.0*M_PI);
  for(int Mu=0; Mu<3; Mu++)
   PAbs1 += PPF*real( Omega*Omega*MMuNu[Mu][Mu]);
  for(int a=0; a<3; a++)
   for(int Mu=0; Mu<3; Mu++)
    PAbs2 += PPF*Omega*Omega*real( conj(p[a][Mu])*p[a][Mu] ); 

  fprintf(f,"%e ",PAbs1);
  fprintf(f,"%e %e %e ",F1[0],F1[1],F1[2]);
  fprintf(f,"%e %e %e ",F2[0],F2[1],F2[2]);
  fprintf(f,"%e %e %e ",F3[0],F3[1],F3[2]);
  fprintf(f,"%e %e %e ",Torque[0],Torque[1],Torque[2]);
  fprintf(f,"%e %e %e ",Fmxp[0], Fmxp[1], Fmxp[2]);
  fprintf(f,"%e %e %e ",Tpxp[0], Tpxp[1], Tpxp[2]);
  fprintf(f,"%e ",PAbs2);

}
