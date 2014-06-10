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

cdouble ExpRel(cdouble x, int n);

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
  zI[0] = ( FA[0]*(Y[0][0]*FB[0] + Y[0][1]*FB[1] + Y[0][2]*FB[2])
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
  int nnz=0;

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

     Entries[nnz] = Result;
     Indices[nnz] = nfB;
     nnz++;
   };

  if (FA->iMTet==-1) return nnz; 

  /*--------------------------------------------------------------*/
  /* handle interactions between face #nfA and face #nfB, where   */
  /* nfB runs over the 3 faces of the negative tetrahedron of nfA */
  /* that differ from nfA                                         */
  /*--------------------------------------------------------------*/
  nt            = FA->iMTet;
  T             = V->Tets[nt];
  Data->QA      = V->Vertices + 3*(FA->iQM);
  Data->PreFacA = -1.0*FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB==nfA) continue;
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

     Entries[nnz] = Result;
     Indices[nnz] = nfB;
     nnz++;

   }; // for(int ifB=0; ifB<4; ifB++)

  return nnz;

} // routine GetVInverse

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GSIData
 { 
   cdouble k;
   double *QA;
   double *QB;

 } GSIData;

void GSurfaceIntegrand(double *xA, double *bA, double DivbA, double *nHatA,
                       double *xB, double *bB, double DivbB, double *nHatB,
                       void *UserData, double *I)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GSIData *Data = (GSIData *) UserData;
  cdouble k     = Data->k;
  double *QA    = Data->QA;
  double *QB    = Data->QB;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double R[3];
  R[0] = (xA[0]-xB[0]);
  R[1] = (xA[1]-xB[1]);
  R[2] = (xA[2]-xB[2]);
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
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
  cdouble Xi3 = Xi*Xi*Xi;
  cdouble w   = ExpRel(Xi, 3);
  cdouble h   = ExpRel(Xi, 2) / Xi;
  cdouble p   = h/Xi - w/Xi3 - 1.0/3.0;

  cdouble T1 = (bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2]) * h;

  cdouble T2 = -DivbA*DivbB*(9.0*h + w + k2*DQdotR*p) / (9.0*k2);

  cdouble *zI = (cdouble *)I;

  zI[0] = NdotN * (T1 + T2) / (-4.0*M_PI*II*k);

}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions using the surface-integral approach of Bleszynski */
/* et. al.                                                     */
/***************************************************************/
cdouble GetGMatrixElement_SI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega,
                             int NumPts=0)
{
  int fDim=2;

  GSIData MyData, *Data = &MyData;

  Data->k = Omega;
  
  SWGFace *FA = VA->Faces[nfA];
  SWGFace *FB = VB->Faces[nfB];

  // 64 surface-surface integrals
  cdouble RetVal=0.0;
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
          double PResult[2], PError[2];
          FaceFaceInt(VA, ntA, nfP, nfBFA, ASign,
                      VB, ntB, nfQ, nfBFB, BSign,
                      GSurfaceIntegrand, (void *)Data, fDim,
                      PResult, PError, NumPts, 10000, 1.0e-8);
          RetVal += cdouble( PResult[0], PResult[1] );
        };
    };

  return RetVal;
}

/***************************************************************/
/* get the dipole and quadupole moments of the current         */
/* distribution described by a single SWG basis function.      */
/***************************************************************/
void GetDQMoments(SWGVolume *O, int nf, double J[3], double Q[3][3],
                  int NeedQ)
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGMatrixElement(SWGVolume *VA, int nfA,
                          SWGVolume *VB, int nfB,
                          cdouble Omega)
{
  SWGFace *FA = VA->Faces[nfA];
  SWGFace *FB = VB->Faces[nfB];

  double rRel = VecDistance(FA->Centroid, FB->Centroid)
                 / fmax(FA->Radius, FB->Radius);

  if ( rRel > 10.0 )
   return GetGMatrixElement_DA(VA, nfA, VB, nfB, Omega, 1);
  else if ( rRel > 2.0 )
   return GetGMatrixElement_DA(VA, nfA, VB, nfB, Omega, 2);
  else 
   return GetGMatrixElement_SI(VA, nfA, VB, nfB, Omega, 9);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AssembleVIEMatrixBlock(int noa, int nob, cdouble Omega,
                                             HMatrix *M, int RowOffset, int ColOffset)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGVolume *OA = Objects[noa];
  SWGVolume *OB = Objects[nob];
  int NFA = OA->NumInteriorFaces;
  int NFB = OB->NumInteriorFaces;
  int SameObject = (noa==nob) ? 1 : 0;

  Log("Assembling U(%i,%i)",noa,nob);
#ifndef USE_OPENMP
  Log(" no multithreading...");
#else
  int NumThreads=GetNumThreads();
  Log(" OpenMP multithreading (%i threads...)",NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nfa=0; nfa<NFA; nfa++)
   for(int nfb=SameObject*nfa; nfb<NFB; nfb++)
    { if (nfb==SameObject*nfa) LogPercent(nfa, NFA);
      M->SetEntry(RowOffset + nfa, ColOffset + nfb,
                  GetGMatrixElement(OA, nfa, OB, nfb, Omega) );
    };

  /***************************************************************/
  /**************************************************************/
  /***************************************************************/
  if (noa==nob)
   { 
     Log("Adding TInv(%i)",noa);
#ifdef USE_OPENMP
       int NumThreads=GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
       for(int nfa=0; nfa<OA->NumInteriorFaces; nfa++)
        { int Indices[7];
          cdouble Entries[7];
          int NNZ = GetVInverseElements(OA, nfa, Omega, Indices, Entries);
          for(int nnz=0; nnz<NNZ; nnz++)
           { int nfb = Indices[nnz];
             M->AddEntry(RowOffset + nfa, ColOffset + nfb, Entries[nnz]);
           };
        };
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AssembleVIEMatrix(cdouble Omega, HMatrix *M)
{
  for(int noa=0; noa<NumObjects; noa++)
   for(int nob=noa; nob<NumObjects; nob++)
    AssembleVIEMatrixBlock(noa, nob, Omega, M, 
                           BFIndexOffset[noa], BFIndexOffset[nob]);

  for(int nr=1; nr<TotalBFs; nr++)
   for(int nc=0; nc<nr; nc++)
    M->SetEntry(nr, nc, conj(M->GetEntry(nc,nr)));

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AllocateVIEMatrix()
{
  return new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX);
}

} // namespace buff
