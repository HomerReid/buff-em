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
  for(int Mu=0; Mu<3; Mu++)
   { vTemp->Copy(vMu[Mu]);
     SInverse->Apply(vTemp, vMu[Mu]);
   };

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { vTemp->Copy(vMuNu[Mu][Nu]);
      SInverse->Apply(vTemp, vMuNu[Mu][Nu]);
    };

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
                  HMatrix *Rytov)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble MMuNu[3][3], MMuNuRho[3][3][3];
  GetNEQMoments(BNEQD, no, Rytov, MMuNu, MMuNuRho);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=vfopen("%s.NEQMoments","a",BNEQD->FileBase);
  fprintf(f,"**MMuNu(Omega=%e): ",Omega);
  fprintf(f," {%+.3e,%+.3e}  {%+.3e,%+.3e }  {%+.3e,%+.3e }\n",
          real(MMuNu[0][0]),imag(MMuNu[0][0]),
          real(MMuNu[0][1]),imag(MMuNu[0][1]),
          real(MMuNu[0][2]),imag(MMuNu[0][2]));
  fprintf(f," {%+.3e,%+.3e}  {%+.3e,%+.3e }  {%+.3e,%+.3e }\n",
          real(MMuNu[1][0]),imag(MMuNu[1][0]),
          real(MMuNu[1][1]),imag(MMuNu[1][1]),
          real(MMuNu[1][2]),imag(MMuNu[1][2]));
  fprintf(f," {%+.3e,%+.3e}  {%+.3e,%+.3e }  {%+.3e,%+.3e }\n",
          real(MMuNu[2][0]),imag(MMuNu[2][0]),
          real(MMuNu[2][1]),imag(MMuNu[2][1]),
          real(MMuNu[2][2]),imag(MMuNu[2][2]));
  fprintf(f,"\n");
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double F1[3], F2[3], F3[3], Torque[3];
  double FPF = TENTHIRDS*ZVAC*Omega*Omega*Omega/(120.0*M_PI);
  for(int Mu=0; Mu<3; Mu++)
   F1[Mu]=F2[Mu]=F3[Mu]=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { 
      F1[Mu] +=     FPF*imag( MMuNuRho[Nu][Mu][Nu] );
      F2[Mu] +=     FPF*imag( MMuNuRho[Mu][Nu][Nu] );
      F3[Mu] -= 4.0*FPF*imag( MMuNuRho[Nu][Nu][Mu] );
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
  f=vfopen("%s.MomentPFT","r",BNEQD->FileBase);
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
  fprintf(f,"\n");
  fclose(f);

}
