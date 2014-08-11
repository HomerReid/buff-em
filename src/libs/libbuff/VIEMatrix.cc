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
/***************************************************************/
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

void VInvIntegrand(double *x, double *b, double DivB, 
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

  cdouble *zI = (cdouble *)I;
  zI[0] = -( FA[0]*(Y[0][0]*FB[0] + Y[0][1]*FB[1] + Y[0][2]*FB[2])
            +FA[1]*(Y[1][0]*FB[0] + Y[1][1]*FB[1] + Y[1][2]*FB[2])
            +FA[2]*(Y[2][0]*FB[0] + Y[2][1]*FB[1] + Y[2][2]*FB[2])
           ) / (Omega*Omega);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetVInverseElements(SWGVolume *V, int nfA, cdouble Omega,
                        int Indices[7], cdouble Entries[7])
{
  Indices[0]=0.0;
  Entries[0]=0.0;
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

     cdouble Result, Error;
     TetInt(V, nt, 0, 1.0, VInvIntegrand, (void *)Data,
            2, (double *)&Result, (double *)&Error, 33, 0, 1.0e-4);

     if (nfB==nfA)
      { 
        Entries[0] += Result;
      }
     else
      { Entries[NNZ] = Result;
        Indices[NNZ] = nfB;
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

     cdouble Result, Error;
     TetInt(V, nt, 0, 1.0, VInvIntegrand, (void *)Data,
            2, (double *)&Result, (double *)&Error, 33, 0, 0);

     if (nfB==nfA)
      { 
        Entries[0] += Result;
      }
     else
      { Entries[NNZ] = Result;
        Indices[NNZ] = nfB;
        NNZ++;
      };

   }; // for(int ifB=0; ifB<4; ifB++)

  return NNZ;

} // routine GetVInverse

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GVIData
 {
   cdouble k;
   bool NeedGradient;
 } GVIData;
 
void GVIntegrand(double *xA, double *bA, double DivbA,
                 double *xB, double *bB, double DivbB,
                 void *UserData, double *I)
{
  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);

  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  if ( fabs(r) < 1.0e-12 ) 
   { I[0]=I[1]=0.0;
     return;
   };

  GVIData *Data     = (GVIData *)UserData;
  cdouble k         = Data->k;
  bool NeedGradient = Data->NeedGradient;

  double DotProduct = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  cdouble PolyFac = DotProduct - DivbA*DivbB/(k*k);

  cdouble IKR = II*k*r;
  cdouble Phi= exp(IKR)/(4.0*M_PI*r);

  cdouble *zI = (cdouble *)I;
  zI[0] = PolyFac*Phi;

  if (NeedGradient)
   { 
     cdouble Psi = (IKR-1.0) * Phi / (r*r);
     zI[1] = R[0]*PolyFac*Psi;
     zI[2] = R[1]*PolyFac*Psi;
     zI[3] = R[2]*PolyFac*Psi;
   };

}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions using a volume cubature method.                   */
/* If GradG is non-null, then GradG[0..2] get dG/dR_i, i=x,y,z */
/***************************************************************/
cdouble GetGMatrixElement_VI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, int Order=0,
                             cdouble *GradG=0)
{
  GVIData MyData, *Data = &MyData;
  Data->k               = Omega;
  Data->NeedGradient    = (GradG!=0);

  int nFun = GradG ? 8 : 2;

  cdouble Result[4], Error[4];
  BFBFInt(VA, nfA, VB, nfB, GVIntegrand, (void *)Data, nFun, 
          (double *)Result, (double *)Error, Order, 10000, 1.0e-8);

  if (GradG)
   { GradG[0] = Result[1];
     GradG[1] = Result[2];
     GradG[2] = Result[3];
   };

  return Result[0];
}

/***************************************************************/
/* routine to compute the h, w, p functions described in the   */
/* memo                                                        */
/***************************************************************/
// h/=Xi;
// cdouble p = h/Xi - w/Xi3 - 1.0/3.0;
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL  
void hwp(cdouble x, cdouble hwp[3], cdouble *hwpPrime)
{
  if ( abs(x) >= 0.1 )
   { 
     cdouble x2=x*x, x3=x2*x, x4=x3*x;
     cdouble ExpRel1, ExpRel2, ExpRel3;

     ExpRel1 = exp(x) - 1.0;
     ExpRel2 = ExpRel1 - x;
     ExpRel3 = ExpRel2 - 0.5*x2;

     hwp[0] = ExpRel2 / x;
     hwp[1] = ExpRel3;
     hwp[2] = hwp[0]/x  - hwp[1]/x3 - 1.0/3.0;

     if (hwpPrime)
      { hwpPrime[0] = ExpRel1/x - ExpRel2/x2;
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

     cdouble Term2=x*x/2.0, Term=Term2, Sum=0.0;
     double nFact = 2.0;
     for(int n=3; n<100; n++)
      { 
        nFact *= n;
        double dn = n;

        hwp[0] += xn  / ( (dn+1.0) * nFact );
        hwp[1] += xn  / nFact;
        hwp[2] += xn / ( nFact*(dn+1.0)*(dn+3.0) );

        if (hwpPrime)
         { hwpPrime[0] += xn / ( (dn+2.0) * nFact);
           hwpPrime[2] += xn / ( (dn+2.0)*(dn+4.0)*nFact);
         };

        Term*=x/((double)m);
        Sum+=Term;
        if ( norm(Term) < EXPRELTOL2*norm(Sum) )
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
  cdouble hwp[3], hwpPrime[3];
  hwp(Xi, hwp, NeedDerivatives ? hwpPrime : 0 );

  cdouble h=hwp[0];
  cdouble w=hwp[1];
  cdouble p=hwp[2];

  cdouble k2  = k*k;
  cdouble Xi  = II*k*r;

  double DotProd = (bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2]);
  double PreFac  = -DivbA*DivbB / (9.0*k2);
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
        T1 = II*k*R[Mu]*DotProd / r; 
        T2 = II*k*R[Mu]*PreFac*(9.0*hP + wP + k2*DQdotR*pP)
              +k2*PreFac*(QA[Mu]-QB[Mu])*p;

        zI[1+Mu] = NdotN * (T1 + T2) / (-4.0*M_PI*II*k);

        zI[3+Mu] = 0.0;
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
                             cdouble Omega, int Order=0,
                             cdouble *dG=0)
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
                      (double *)Result, (double *)Error, Order, 10000, 1.0e-8);
          RetVal += Result[0];
          if (dG)
           for(int n=0; n<6; n++) dG[n]+=Result[n+1];
        };
    };

  return RetVal;
}

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
   return GetGMatrixElement_SI(VA, nfA, VB, nfB, Omega, 20, dG);
  else if (ncv==1)
   return GetGMatrixElement_SI(VA, nfA, VB, nfB, Omega, 9,  dG);
  else
   return GetGMatrixElement_VI(VA, nfA, VB, nfB, Omega, 16, dG);
  
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
  cdouble dGBuffer[6];
  cdouble *dG = dGMatrix==0 ? 0 : dGBuffer;

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
      G->SetEntry(Row, Col, GetGMatrixElement(OA, nfa, OB, nfb, Omega, dG) );

      if (dG)
       for(int Mu=0; Mu<6; Mu++)
        if (dGMatrix[Mu]) 
         dGMatrix[Mu]->SetEntry(Row, Col, dG[Mu]);
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleVInvBlock(int no, cdouble Omega,
                                    SMatrix *VInv, SMatrix *ImV,
                                    HMatrix *TInv=0, int Offset=0)
{
   Log("Adding TInv(%i)",no);
   SWGObject *O = Objects[no];
   int Offset = BFIndexOffset[no];

//BeginAssembly(int est_nnz = 0 /* default: allocate on fly */);

#ifdef USE_OPENMP
   int NumThreads=GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
   for(int nr=0; nf<O->NumInteriorFaces; nr++)
    { int nc[7];
      cdouble VInvEntries[7];
      double ImvVntries[7];
      int NNZ = GetVInvElements(O, nr, Omega, nc,
                                VInvEntries, ImVEntries);
      for(int nnz=0; nnz<NNZ; nnz++)
        { 
          if (VInv) 
           VInv->SetEntry(nr, nc[nnz], VInvEntries[nnz]);
          if (ImV) 
           ImV->SetEntry(nr, nc[nnz], ImVEntries[nnz]);
          if (TInv)
           TInv->AddEntry(Offset + nr, Offset + nc[nnz], VInvEntries[nnz]);
        };
   };

//if (VInv)
// VInv->EndAssembly

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleVBlock(int noa, int nob, cdouble Omega, 
                                 HMatrix *U, HMatrix **GradU)
{ AssembleVIEMatrixBlock(noa, nob, Omega, U, GradU); }

void SWGGeometry::AssembleTInvBlock(int no, cdouble Omega, HMatrix *TInv)
{ AssembleVIEMatrixBlock(no, no, Omega, TInv); }

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
       AssembleVInvBlock(noa, Omega, 0, 0, M, 
                         BFIndexOffset[noa], BFIndexOffset[noa]);
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
