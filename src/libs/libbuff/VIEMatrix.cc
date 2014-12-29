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
 * VIEMatrix.cc -- libSWG routines for assembling matrices and vectors
 *
 * homer reid   -- 5/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libSGJC.h"
#include "libscuff.h"
#include "libbuff.h"

namespace scuff{

void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

void CalcGC(double R1[3], double R2[3],
            cdouble Omega, cdouble EpsR, cdouble MuR, 
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

}

using namespace scuff;

namespace buff {

#define MAXSTR 1000
#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Invert3x3Matrix(cdouble M[3][3], cdouble W[3][3])
{
  cdouble M11=M[0][0], M12=M[0][1], M13=M[0][2];
  cdouble M21=M[1][0], M22=M[1][1], M23=M[1][2];
  cdouble M31=M[2][0], M32=M[2][1], M33=M[2][2];

  cdouble Det =  M12*M23*M31 + M13*M21*M32 + M11*M22*M33
                -M13*M22*M31 - M11*M23*M32 - M12*M21*M33;

  W[0][0] = (M22*M33 - M23*M32) / Det;
  W[0][1] = (M13*M32 - M12*M33) / Det;
  W[0][2] = (M12*M23 - M13*M22) / Det;

  W[1][0] = (M23*M31 - M21*M33) / Det;
  W[1][1] = (M11*M33 - M13*M31) / Det;
  W[1][2] = (M13*M21 - M11*M23) / Det;

  W[2][0] = (M21*M32 - M22*M31) / Det;
  W[2][1] = (M12*M31 - M11*M32) / Det;
  W[2][2] = (M11*M22 - M12*M21) / Det;

}

/***************************************************************/
/* user data structure and integrand function for **************/
/* GetVInvAndImEpsEntries                                      */
/***************************************************************/
typedef struct VIIData
 { 
   double *QA;
   double PreFacA;
   double *QB;
   double PreFacB;
   cdouble Omega;
   IHAIMatProp *MP;

 } VIIData;

void VIntegrand(double *x, double *b, double DivB, 
                void *UserData, double *I)
{
  VIIData *Data   = (VIIData *)UserData;
  double *QA      = Data->QA;
  double PreFacA  = Data->PreFacA;
  double *QB      = Data->QB;
  double PreFacB  = Data->PreFacB;
  cdouble Omega   = Data->Omega;
  IHAIMatProp *MP = Data->MP;

  cdouble Eps[3][3], Y[3][3];
  MP->GetEps( Omega, x, Eps );

  Eps[0][0] -= 1.0;
  Eps[1][1] -= 1.0;
  Eps[2][2] -= 1.0;
  Invert3x3Matrix(Eps, Y);

  double FA[3], FB[3];
  FA[0] = PreFacA * (x[0] - QA[0]);
  FA[1] = PreFacA * (x[1] - QA[1]);
  FA[2] = PreFacA * (x[2] - QA[2]);
  FB[0] = PreFacB * (x[0] - QB[0]);
  FB[1] = PreFacB * (x[1] - QB[1]);
  FB[2] = PreFacB * (x[2] - QB[2]);

  cdouble VInvElement=0.0;
  double ImEpsElement=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { VInvElement  += FA[Mu]*Y[Mu][Nu]*FB[Nu];
      ImEpsElement += FA[Mu]*imag(Eps[Mu][Nu])*FB[Nu];
    };
  VInvElement *= -1.0/ (Omega*Omega);
 
  I[0] = real(VInvElement);
  I[1] = imag(VInvElement);
  I[2] = ImEpsElement;
  
}

/***************************************************************/
/* For a given SWG basis function f_a, this routine computes   */
/* the matrix elements of the VInverse and Im Eps operators    */
/* (where  VInverse=Eps^{-1} / k^2) between f_a and all basis  */
/* functions f_b that have nonzero overlap with f_a (there are */
/* a maximum of 7 such BFs.) The indices of the f_b functions  */
/* are returned in Indices.                                    */
/***************************************************************/
int GetVInvAndImEpsEntries(SWGVolume *V, int nfA, cdouble Omega, int Indices[7],
                           cdouble VInvEntries[7], double ImEpsEntries[7])
{
  Indices[0]=0.0;
  VInvEntries[0]=0.0;
  ImEpsEntries[0]=0.0;
  int NNZ=1;

  SWGFace *FA = V->Faces[nfA];
  struct VIIData MyVIIData, *Data=&MyVIIData;
  Data->Omega = Omega;
  Data->MP    = V->MP;

  /*--------------------------------------------------------------*/
  /* handle interactions between face #nfA and face #nfB, where   */
  /* nfB runs over the 4 faces of the positive tetrahedron of nfA */
  /*--------------------------------------------------------------*/
  int nt        = FA->iPTet;
  SWGTet *T     = V->Tets[nt];
  Data->QA      = V->Vertices + 3*(FA->iQP);
  Data->PreFacA = FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB >= V->NumInteriorFaces) continue;

     SWGFace *FB = V->Faces[nfB];
     if ( FB->iPTet == nt )
      { Data->QB      = V->Vertices + 3*(FB->iQP);
        Data->PreFacB = FB->Area / (3.0*T->Volume);
      }
     else
      { Data->QB      = V->Vertices + 3*(FB->iQM);
        Data->PreFacB = -1.0*FB->Area / (3.0*T->Volume);
      };

     double I[3], E[3];
     TetInt(V, nt, 0, 1.0, VIntegrand, (void *)Data,
            3, I, E, 33, 0, 1.0e-4);

     if (nfB==nfA)
      { 
        VInvEntries[0]  += cdouble(I[0], I[1]);
        ImEpsEntries[0] += I[2];
      }
     else
      { VInvEntries[NNZ]  = cdouble(I[0], I[1]);
        ImEpsEntries[NNZ] = I[2];
        Indices[NNZ]      = nfB;
        NNZ++;
      };

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  nt            = FA->iMTet;
  T             = V->Tets[nt];
  Data->QA      = V->Vertices + 3*(FA->iQM);
  Data->PreFacA = -1.0*FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB >= V->NumInteriorFaces) continue;

     SWGFace *FB = V->Faces[nfB];
     if ( FB->iPTet == nt )
      { Data->QB      = V->Vertices + 3*(FB->iQP);
        Data->PreFacB = FB->Area / (3.0*T->Volume);
      }
     else
      { Data->QB      = V->Vertices + 3*(FB->iQM);
        Data->PreFacB = -1.0*FB->Area / (3.0*T->Volume);
      };

     double I[3], E[3];
     TetInt(V, nt, 0, 1.0, VIntegrand, (void *)Data,
            3, I, E, 33, 0, 0);

     if (nfB==nfA)
      { 
        VInvEntries[0]   += cdouble( I[0], I[1] );
        ImEpsEntries[0]  += I[2];
      }
     else
      { VInvEntries[NNZ]  = cdouble( I[0], I[1] );
        ImEpsEntries[NNZ] = I[2];
        Indices[NNZ] = nfB;
        NNZ++;
      };

   }; // for(int ifB=0; ifB<4; ifB++)

  return NNZ;

} // routine GetVInvAndImEpsEntries

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GVIData
 {
   cdouble k;
   double x0[3]; // origin for torque calculations
   bool NeedDerivatives;
 } GVIData;
 
void GVIntegrand(double *xA, double *bA, double DivbA,
                 double *xB, double *bB, double DivbB,
                 void *UserData, double *I)
{
  GVIData *Data        = (GVIData *)UserData;
  bool NeedDerivatives = Data->NeedDerivatives;
  cdouble k            = Data->k;

  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);

  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  if ( fabs(r) < 1.0e-12 ) 
   { if (Data->NeedDerivatives)
      memset(I,0,14*sizeof(double)); 
     else
      I[0]=I[1]=0.0;
     return;
   };

  double DotProduct = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  cdouble PolyFac = DotProduct - DivbA*DivbB/(k*k);

  cdouble IKR = II*k*r;
  cdouble Phi= exp(IKR)/(4.0*M_PI*r);

  cdouble *zI = (cdouble *)I;
  zI[0] = PolyFac*Phi;

  if (NeedDerivatives)
   { 
     cdouble Psi = (IKR-1.0) * Phi / (r*r);
     zI[1] = R[0]*PolyFac*Psi;
     zI[2] = R[1]*PolyFac*Psi;
     zI[3] = R[2]*PolyFac*Psi;

     double xAmx0[3];
     double *x0 = Data->x0;
     xAmx0[0] = xA[0] - x0[0];
     xAmx0[1] = xA[1] - x0[1];
     xAmx0[2] = xA[2] - x0[2];
     for(int Mu=0; Mu<3; Mu++)
      { int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
        zI[4 + Mu] 
         = (xAmx0[MP1]*bB[MP2]-xAmx0[MP2]*bB[MP1])*Phi*DivbA/3.0
          +(xAmx0[MP1]* R[MP2]-xAmx0[MP2]* R[MP1])*PolyFac*Psi;
      };

   };

}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions using a volume cubature method.                   */
/* If dG is non-null, then on return we have                   */
/*  dG[0..2] = dG/dR_i,    i=x,y,z,                            */
/*  dG[3..5] = dG/dTheta_i                                     */
/* where dG/dTheta_i is the angular derivative with respect to */
/* rotations about the ith axis                                */
/***************************************************************/
cdouble GetGMatrixElement_VI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, cdouble *dG=0,
                             int Order=0, int MaxEvals=10000)
{
  GVIData MyData, *Data = &MyData;
  Data->k               = Omega;
  int nFun;

  if (dG==0)
   { Data->NeedDerivatives = false;
     nFun = 2;
   }
  else
   { Data->NeedDerivatives = true;
     nFun = 14;
     Data->x0[0] = Data->x0[1] = Data->x0[2] = 0.0;
     if (VA->OTGT) VA->OTGT->Apply( Data->x0 );
     if (VA->GT) VA->GT->Apply( Data->x0 );
   };

  cdouble Result[7], Error[7];
  BFBFInt(VA, nfA, VB, nfB, GVIntegrand, (void *)Data, nFun,
          (double *)Result, (double *)Error, Order, MaxEvals, 1.0e-8);

  if (dG) memcpy(dG, Result+1, 6*sizeof(cdouble));

  return Result[0];
}

/***************************************************************/
/* routine to compute the h, w, p functions described in the   */
/* memo                                                        */
/***************************************************************/
#define EXPRELTOL  1.0e-8
void Gethwp(cdouble x, cdouble hwp[3], cdouble *hwpPrime)
{
  if ( abs(x) >= 0.1 )
   { 
     cdouble x2=x*x;
     cdouble x3=x2*x;
     cdouble ExpRel1 = exp(x) - 1.0;
     cdouble ExpRel2 = ExpRel1 - x;
     cdouble ExpRel3 = ExpRel2 - 0.5*x2;

     hwp[0] = ExpRel2 / x;
     hwp[1] = ExpRel3;
     hwp[2] = hwp[0]/x  - hwp[1]/x3 - 1.0/3.0;

     if (hwpPrime)
      { cdouble x4=x3*x;
        hwpPrime[0] = ExpRel1/x - ExpRel2/x2;
        hwpPrime[1] = ExpRel2;
        hwpPrime[2] = hwpPrime[0]/x - hwp[0]/x2 + 3.0*hwp[1]/x4 - hwpPrime[1]/x3;
      };
   }
  else
   { 
     cdouble x2=x*x;

     hwp[0] = x/2.0 + x2/6.0;
     hwp[1] = 0.0;
     hwp[2] = x/8.0 + x2/30.0;

     if (hwpPrime)
      { hwpPrime[0] = 1.0/2.0 + x/3.0 + x2/8.0;
        hwpPrime[2] = 1.0/8.0 + x/15.0 + x2/48.0;
      };

     // in this loop, Term = x^n / n!
     cdouble Term = x2 / 2.0;
     for(double n=3.0; n<100.1; n+=1.0)
      { 
        Term *= x / n;

        hwp[0] += Term / (n+1.0);
        hwp[1] += Term;
        hwp[2] += Term / ( (n+3.0)*(n+1.0) );

        if (hwpPrime)
         { hwpPrime[0] += Term / (n+2.0);
           hwpPrime[2] += Term / ( (n+4.0)*(n+2.0) );
         };

        if ( abs(Term) < EXPRELTOL*abs(hwp[1]) )
         break;

      };

     if (hwpPrime)
      hwpPrime[1] = hwp[1] + x2/2.0;

   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GSIData
 { 
   cdouble k;
   double *QA;
   double *QB;
   bool NeedDerivatives;

 } GSIData;

void GSurfaceIntegrand(double *xA, double *bA, double DivbA, double *nHatA,
                       double *xB, double *bB, double DivbB, double *nHatB,
                       void *UserData, double *I)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GSIData *Data        = (GSIData *) UserData;
  cdouble k            = Data->k;
  double *QA           = Data->QA;
  double *QB           = Data->QB;
  bool NeedDerivatives = Data->NeedDerivatives;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double R[3];
  R[0] = (xA[0]-xB[0]);
  R[1] = (xA[1]-xB[1]);
  R[2] = (xA[2]-xB[2]);
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  if ( fabs(r2) < 1.0e-15 )
   { int fdim = NeedDerivatives ? 14 : 2;
     memset(I, 0, fdim*sizeof(double));
     return;
   };
  double r = sqrt(r2);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double NdotN =   nHatA[0]*nHatB[0] 
                 + nHatA[1]*nHatB[1] 
                 + nHatA[2]*nHatB[2];

  double DQdotR =  (QA[0]-QB[0])*R[0] 
                  +(QA[1]-QB[1])*R[1] 
                  +(QA[2]-QB[2])*R[2];

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble k2  = k*k;
  cdouble Xi  = II*k*r;

  cdouble hwp[3], hwpPrime[3];
  Gethwp(Xi, hwp, NeedDerivatives ? hwpPrime : 0 );

  cdouble h=hwp[0];
  cdouble w=hwp[1];
  cdouble p=hwp[2];

  double DotProd = (bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2]);
  cdouble PreFac  = -DivbA*DivbB / (9.0*k2);
  cdouble T1 = DotProd * h;
  cdouble T2 = PreFac * (9.0*h + w + k2*DQdotR*p);

  cdouble *zI = (cdouble *)I;
  zI[0] = NdotN * (T1 + T2) / (-4.0*M_PI*II*k);

  if (NeedDerivatives)
   { 
     cdouble hP=hwpPrime[0];
     cdouble wP=hwpPrime[1];
     cdouble pP=hwpPrime[2];

     for(int Mu=0; Mu<3; Mu++)
      { 
        cdouble Factor = II*k*R[Mu]/r;

        T1 = Factor*DotProd*hP;
        T2 = Factor*PreFac*(9.0*hP + wP + k2*DQdotR*pP)
              +k2*PreFac*(QA[Mu]-QB[Mu])*p;
        zI[1+Mu] = NdotN * (T1 + T2) / (-4.0*M_PI*II*k);

        zI[4+Mu] = 0.0;
      };
   };

}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions using the surface-integral approach of Bleszynski */
/* et. al.                                                     */
/***************************************************************/
cdouble GetGMatrixElement_SI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, cdouble *dG=0,
                             int Order=0, int MaxEvals=10000)
{

  GSIData MyData, *Data = &MyData;

  Data->k = Omega;
  Data->NeedDerivatives = (dG!=0);

  int fDim = (dG==0) ? 2:14; 
  
  SWGFace *FA = VA->Faces[nfA];
  SWGFace *FB = VB->Faces[nfB];

  // 64 surface-surface integrals
  cdouble RetVal=0.0;
  if (dG) memset(dG, 0, 6*sizeof(cdouble));
  for(int ASign=+1; ASign>=-1; ASign-=2)
   for(int BSign=+1; BSign>=-1; BSign-=2)
    {
      int ntA    = (ASign==1) ? FA->iPTet  : FA->iMTet;
      int nfBFA  = (ASign==1) ? FA->PIndex : FA->MIndex;
      SWGTet *TA = VA->Tets[ ntA ];
      Data->QA   = VA->Vertices + 3*TA->VI[nfBFA];

      int ntB    = (BSign==1) ? FB->iPTet  : FB->iMTet;
      int nfBFB  = (BSign==1) ? FB->PIndex : FB->MIndex;
      SWGTet *TB = VB->Tets[ ntB ];
      Data->QB   = VB->Vertices + 3*TB->VI[nfBFB];

      for(int nfP=0; nfP<4; nfP++)
       for(int nfQ=0; nfQ<4; nfQ++)
        { 
          cdouble Result[7], Error[7];
          FaceFaceInt(VA, ntA, nfP, nfBFA, ASign,
                      VB, ntB, nfQ, nfBFB, BSign,
                      GSurfaceIntegrand, (void *)Data, fDim,
                      (double *)Result, (double *)Error, 
                      Order, MaxEvals, 1.0e-8);
          RetVal += Result[0];

          if (dG) for(int n=0; n<6; n++) dG[n]+=Result[n+1];
        };
    };

  return RetVal;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGMatrixElement(SWGVolume *VA, int nfA,
                          SWGVolume *VB, int nfB,
                          cdouble Omega, cdouble *dG=0)
{
  SWGFace *FA = VA->Faces[nfA];
  SWGFace *FB = VB->Faces[nfB];

  double rRel;
  int ncv = CompareBFs(VA, nfA, VB, nfB, &rRel);

  if ( ncv>=2 )
   return GetGMatrixElement_SI(VA, nfA, VB, nfB, Omega, dG, 20);
  else if (ncv==1)
   return GetGMatrixElement_SI(VA, nfA, VB, nfB, Omega, dG, 9);
  else
   return GetGMatrixElement_VI(VA, nfA, VB, nfB, Omega, dG, 16);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleGBlock(int noa, int nob, cdouble Omega,
                                 HMatrix *G, HMatrix **dGMatrix,
                                 int RowOffset, int ColOffset)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGVolume *OA = Objects[noa];
  SWGVolume *OB = Objects[nob];
  int NFA = OA->NumInteriorFaces;
  int NFB = OB->NumInteriorFaces;
  int SameObject = (noa==nob) ? 1 : 0;

  Log("Assembling G(%i,%i)",noa,nob);
#ifndef USE_OPENMP
  Log(" no multithreading...");
#else
  int NumThreads=GetNumThreads();
  Log(" OpenMP multithreading (%i threads...)",NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nfa=0; nfa<NFA; nfa++)
   for(int nfb=SameObject*nfa; nfb<NFB; nfb++)
    { 
      if (nfb==SameObject*nfa) 
       LogPercent(nfa, NFA);

      int Row=RowOffset + nfa;
      int Col=ColOffset + nfb;
      if (dGMatrix==0)
       {
         G->SetEntry(Row, Col, GetGMatrixElement(OA, nfa, OB, nfb, Omega, 0) );
       }
      else 
       { cdouble dG[6];
         G->SetEntry(Row, Col, GetGMatrixElement(OA, nfa, OB, nfb, Omega, dG) );
         for(int Mu=0; Mu<6; Mu++)
          if (dGMatrix[Mu]) 
           dGMatrix[Mu]->SetEntry(Row, Col, dG[Mu]);
       };
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleVInvBlock(int no, cdouble Omega,
                                    SMatrix *VInv, SMatrix *ImEps,
                                    HMatrix *TInv, int Offset)
{
   Log("Adding VInv(%i)",no);
   SWGVolume *O = Objects[no];

#ifdef USE_OPENMP
   int NumThreads=GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
   for(int nr=0; nr<O->NumInteriorFaces; nr++)
    { int nc[7];
      cdouble VInvEntries[7];
      double ImEpsEntries[7];
      int NNZ = GetVInvAndImEpsEntries(O, nr, Omega, nc, VInvEntries, ImEpsEntries);
      for(int nnz=0; nnz<NNZ; nnz++)
        { 
          if (VInv) 
           VInv->SetEntry(nr, nc[nnz], VInvEntries[nnz]);
          if (ImEps) 
           ImEps->SetEntry(nr, nc[nnz], ImEpsEntries[nnz]);
          if (TInv)
           TInv->AddEntry(Offset + nr, Offset + nc[nnz], VInvEntries[nnz]);
        };
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AssembleVIEMatrix(cdouble Omega, HMatrix *M)
{
  if ( M && ( (M->NR!=TotalBFs) || (M->NR != M->NC) ) )
   { Warn("wrong-size M-matrix passed to AssembleVIEMatrix (reallocating...)");
     delete M;
     M=0;
   };
  if (!M)
   M = new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX);

  for(int noa=0; noa<NumObjects; noa++)
   for(int nob=noa; nob<NumObjects; nob++)
    { 
      AssembleGBlock(noa, nob, Omega, M, 0,
                     BFIndexOffset[noa], BFIndexOffset[nob]);

      if (nob==noa)
       AssembleVInvBlock(noa, Omega, 0, 0, M, BFIndexOffset[noa]);
    };

  for(int nr=1; nr<TotalBFs; nr++)
   for(int nc=0; nc<nr; nc++)
    M->SetEntry(nr, nc, M->GetEntry(nc,nr));

  return M;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AllocateVIEMatrix(bool PureImagFreq)
{
  if (PureImagFreq)
   return new HMatrix(TotalBFs, TotalBFs, LHM_REAL);
  else
   return new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX);
}

} // namespace buff
