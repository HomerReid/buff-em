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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

using namespace scuff;

namespace buff {

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
#define BOLTZMANNK 4.36763e-4
double GetThetaFactor(double Omega, double T)
{ 
  if (T==0.0)
   return 0.0;
  return Omega / ( exp( Omega/(BOLTZMANNK*T) ) - 1.0 );
}

double GetAverageTemperature(SWGVolume *O, SVTensor *TemperatureSVT)
{
  if (TemperatureSVT==0)
   return 0.0;

  double TAvg=0.0, TotalVolume=0.0;
  for(int nt=0; nt<O->NumTets; nt++)
   { 
     double *X0   = O->Tets[nt]->Centroid;
     double V     = O->Tets[nt]->Volume;
     double T     = TemperatureSVT->EvaluateD(0.0, X0);
     TotalVolume += V;
     TAvg        += V*T;
   };

  return TAvg/TotalVolume;

}


/***************************************************************/
/* user data structure and integrand function for              */
/* GetOverlaps                                                 */
/***************************************************************/
typedef struct GOData
 { 
   double *QA;
   double PreFacA;
   double *QB;
   double PreFacB;
   cdouble Omega;
   SVTensor *EpsSVT;
   SVTensor *TemperatureSVT;
   double ThetaEnvironment;
   double DeltaThetaHat;
 } GOData;

void GetOverlapIntegrand(double *x, double *b, double DivB,
                         void *UserData, double *I)
{
  (void) DivB;
  (void) b;
 
  GOData *Data             = (GOData *)UserData;
  double *QA               = Data->QA;
  double PreFacA           = Data->PreFacA;
  double *QB               = Data->QB;
  double PreFacB           = Data->PreFacB;
  double Omega             = real(Data->Omega);
  SVTensor *EpsSVT         = Data->EpsSVT;
  SVTensor *TemperatureSVT = Data->TemperatureSVT;
  double ThetaEnvironment  = Data->ThetaEnvironment;
  double DeltaThetaHat     = Data->DeltaThetaHat;

  cdouble EpsM1[3][3], InvEpsM1[3][3];
  EpsSVT->Evaluate( Omega, x, EpsM1 );
  EpsM1[0][0] -= 1.0;
  EpsM1[1][1] -= 1.0;
  EpsM1[2][2] -= 1.0;
  Invert3x3Matrix(EpsM1, InvEpsM1);

  double FA[3], FB[3];
  FA[0] = PreFacA * (x[0] - QA[0]);
  FA[1] = PreFacA * (x[1] - QA[1]);
  FA[2] = PreFacA * (x[2] - QA[2]);
  FB[0] = PreFacB * (x[0] - QB[0]);
  FB[1] = PreFacB * (x[1] - QB[1]);
  FB[2] = PreFacB * (x[2] - QB[2]);

  double DeltaTheta=1.0;
  if (TemperatureSVT)
   { 
     double T = TemperatureSVT->EvaluateD(0,x);
     DeltaTheta = GetThetaFactor( Omega, T ) - ThetaEnvironment;
   };

  cdouble V=0.0, VInv=0.0;
  double Rytov=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { V     += FA[Mu]*EpsM1[Mu][Nu]*FB[Nu];
      VInv  += FA[Mu]*InvEpsM1[Mu][Nu]*FB[Nu];
      Rytov += FA[Mu]*imag(EpsM1[Mu][Nu])*FB[Nu];
    };
  V     *= -1.0*Omega*Omega;
  VInv  *= -1.0/(Omega*Omega);
  Rytov *= 2.0*Omega*DeltaTheta/(M_PI*ZVAC*DeltaThetaHat);
 
  I[0] = real(V);
  I[1] = imag(V);
  I[2] = real(VInv);
  I[3] = imag(VInv);
  I[4] = Rytov; 
  
}

/***************************************************************/
/* For a given SWG basis function f_a, this routine computes   */
/* the matrix elements of various operators between f_a and    */
/* all basis functions f_b that have nonzero overlap with f_a  */
/* (there are at most MAXOVERLAP=7 such BFs.)                  */
/* The indices of the f_b functions are returned in Indices.   */
/*                                                             */
/*        V = k^2 * (1 - Eps)k^2                               */
/* VInverse = V^{-1}                                           */
/*    Rytov = (2k / Pi*ZVAC) * (RelDelTheta(T)) * Im Eps       */
/*                                                             */
/* where RelDelTheta is defined as follows:                    */
/*  DelTheta(x)    = Theta( T(x) ) - Theta( TEnvironment )     */
/*  DelThetaHat(x) = Theta( TAvg ) - Theta( TEnvironment )     */
/*  RelDelTheta    = DelTheta(x) / DeltaThetaHat               */
/*                                                             */
/* If TemperatureSVT==NULL then the factor Theta(T)-Theta(TEnv)*/
/* is replaced by 1.                                           */
/***************************************************************/
int GetOverlaps(SWGVolume *O, int nfA, cdouble Omega,
                SVTensor *TemperatureSVT,
                double TAvg, double TEnvironment,
                int Indices[MAXOVERLAP],
                cdouble VEntries[MAXOVERLAP],
                cdouble VInvEntries[MAXOVERLAP],
                double RytovEntries[MAXOVERLAP])
{
  Indices[0]=nfA;
  VEntries[0]=0.0;
  VInvEntries[0]=0.0;
  RytovEntries[0]=0.0;
  int NNZ=1;

  double ThetaAvg         = GetThetaFactor( real(Omega), TAvg);
  double ThetaEnvironment = GetThetaFactor( real(Omega), TEnvironment);

  SWGFace *FA = O->Faces[nfA];
  struct GOData MyGOData, *Data=&MyGOData;
  Data->Omega            = Omega;
  Data->EpsSVT           = O->SVT;
  Data->TemperatureSVT   = TemperatureSVT;
  Data->ThetaEnvironment = ThetaEnvironment;
  Data->DeltaThetaHat    = ThetaAvg - ThetaEnvironment;

  /*--------------------------------------------------------------*/
  /* handle interactions between face #nfA and face #nfB, where   */
  /* nfB runs over the 4 faces of the positive tetrahedron of nfA */
  /*--------------------------------------------------------------*/
  int nt        = FA->iPTet;
  SWGTet *T     = O->Tets[nt];
  Data->QA      = O->Vertices + 3*(FA->iQP);
  Data->PreFacA = FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB >= O->NumInteriorFaces) continue;

     SWGFace *FB = O->Faces[nfB];
     if ( FB->iPTet == nt )
      { Data->QB      = O->Vertices + 3*(FB->iQP);
        Data->PreFacB = FB->Area / (3.0*T->Volume);
      }
     else
      { Data->QB      = O->Vertices + 3*(FB->iQM);
        Data->PreFacB = -1.0*FB->Area / (3.0*T->Volume);
      };

     #define NFUN 5
     int Order=33;
     double RelTol=1.0e-4;
     double I[NFUN], E[NFUN];
     TetInt(O, nt, 0, 1.0, GetOverlapIntegrand, (void *)Data,
            NFUN, I, E, Order, 0, RelTol);

     if (nfB==nfA)
      { VEntries[0]     += cdouble(I[0], I[1]);
        VInvEntries[0]  += cdouble(I[2], I[3]);
        RytovEntries[0] += I[4];
      }
     else
      { VEntries[NNZ]     = cdouble(I[0], I[1]);
        VInvEntries[NNZ]  = cdouble(I[2], I[3]);
        RytovEntries[NNZ] = I[4];
        Indices[NNZ]      = nfB;
        NNZ++;
      };

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  nt            = FA->iMTet;
  T             = O->Tets[nt];
  Data->QA      = O->Vertices + 3*(FA->iQM);
  Data->PreFacA = -1.0*FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB >= O->NumInteriorFaces) continue;

     SWGFace *FB = O->Faces[nfB];
     if ( FB->iPTet == nt )
      { Data->QB      = O->Vertices + 3*(FB->iQP);
        Data->PreFacB = FB->Area / (3.0*T->Volume);
      }
     else
      { Data->QB      = O->Vertices + 3*(FB->iQM);
        Data->PreFacB = -1.0*FB->Area / (3.0*T->Volume);
      };

     int Order=33;
     double RelTol=1.0e-4;
     double I[NFUN], E[NFUN];
     TetInt(O, nt, 0, 1.0, GetOverlapIntegrand, (void *)Data,
            NFUN, I, E, Order, 0, RelTol);

     if (nfB==nfA)
      { 
        VEntries[0]     += cdouble( I[0], I[1] );
        VInvEntries[0]  += cdouble( I[2], I[3] );
        RytovEntries[0] += I[4];
      }
     else
      { VEntries[NNZ]     = cdouble( I[0], I[1] );
        VInvEntries[NNZ]  = cdouble( I[2], I[3] );
        RytovEntries[NNZ] = I[4];
        Indices[NNZ]      = nfB;
        NNZ++;
      };

   }; // for(int ifB=0; ifB<4; ifB++)

  return NNZ;

} // routine GetOverlaps

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleOverlapBlocks(int no, cdouble Omega,
                                        SVTensor *TemperatureSVT,
                                        double TAvg,
                                        double TEnvironment,
                                        SMatrix *V,
                                        SMatrix *VInv,
                                        SMatrix *Rytov,
                                        HMatrix *TInv,
                                        int Offset)
{
   SWGVolume *O = Objects[no];

#ifdef USE_OPENMP
   int NumThreads=GetNumThreads();
   if (V || VInv || Rytov)
    NumThreads=1;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
   for(int nr=0; nr<O->NumInteriorFaces; nr++)
    { int ncList[MAXOVERLAP];
      cdouble VEntries[MAXOVERLAP];
      cdouble VInvEntries[MAXOVERLAP];
      double RytovEntries[MAXOVERLAP];
      int NNZ = GetOverlaps(O, nr, Omega,
                            TemperatureSVT, TAvg, TEnvironment,
                            ncList, VEntries, VInvEntries, RytovEntries);
      for(int nnz=0; nnz<NNZ; nnz++)
        { int nc = ncList[nnz];
          if (V)
           V->SetEntry(nr, nc, VEntries[nnz]);
          if (VInv)
           VInv->SetEntry(nr, nc, VInvEntries[nnz]);
          if (Rytov)
           Rytov->SetEntry(nr, nc, RytovEntries[nnz]);
          if (TInv)
           TInv->AddEntry(Offset+nr, Offset+nc, VInvEntries[nnz]);
        };
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleGBlock(int noa, int nob, cdouble Omega,
                                 HMatrix *G,
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
  Log("AGB Assembling G(%i,%i)",noa,nob);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIBBICache *GCache   = 0;
  if (noa==nob)
   {
     int noMate=Mate[noa];
     int noCache = (noMate==-1) ? noa : noMate;
     if (ObjectGCaches[noCache]==0)
      ObjectGCaches[noCache]=new FIBBICache(OA->MeshFileName);
     GCache = ObjectGCaches[noCache];
     Log("AGB initial cache size %i ",GCache->Size());
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#ifndef USE_OPENMP
  Log("AGB no multithreading...");
#else
  int NumThreads=GetNumThreads();
  Log("AGB OpenMP multithreading (%i threads...)",NumThreads);
#pragma omp parallel for schedule(dynamic,1),		\
                         num_threads(NumThreads)
#endif
  for(int nfa=0; nfa<NFA; nfa++)
   for(int nfb=SameObject*nfa; nfb<NFB; nfb++)
    { 
      if (nfb==SameObject*nfa)
       LogPercent(nfa, NFA);

      int Row=RowOffset + nfa;
      int Col=ColOffset + nfb;
      cdouble GAB=GetGMatrixElement(OA, nfa, OB, nfb, Omega, GCache);
      G->SetEntry(Row, Col, GAB);
      if (SameObject && nfb>nfa)
       G->SetEntry(Col, Row, GAB);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (GCache)
   { Log("AGB final cache size %i ",GCache->Size());
     GCache->Store(OA->MeshFileName);
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
      AssembleGBlock(noa, nob, Omega, M, 
                     BFIndexOffset[noa], BFIndexOffset[nob]);

      if (nob==noa)
       { Log("Adding VInv(%i)",noa);
         AssembleOverlapBlocks(noa, Omega, 0, 0.0, 0.0, 0, 0, 0,
                               M, BFIndexOffset[noa]);
 
       };
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
