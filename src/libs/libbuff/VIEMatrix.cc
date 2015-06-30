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

/***************************************************************/
/* user data structure and integrand function for              */
/* GetVInvAndSigmaEntries                                 */
/***************************************************************/
typedef struct VIIData
 { 
   double *QA;
   double PreFacA;
   double *QB;
   double PreFacB;
   cdouble Omega;
   IHAIMatProp *MP;
   IHAIMatProp *Temperature;

 } VIIData;

void VIntegrand(double *x, double *b, double DivB, 
                void *UserData, double *I)
{
  (void) DivB;
  (void) b;

  VIIData *Data            = (VIIData *)UserData;
  double *QA               = Data->QA;
  double PreFacA           = Data->PreFacA;
  double *QB               = Data->QB;
  double PreFacB           = Data->PreFacB;
  cdouble Omega            = Data->Omega;
  IHAIMatProp *MP          = Data->MP;
  IHAIMatProp *Temperature = Data->Temperature;

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

  double Theta=1.0;
  if (Temperature)
   { 
     cdouble TT[3][3];
     Temperature->GetEps(0,x,TT);
     double T=real(TT[0][0]);
     Theta=GetThetaFactor(real(Omega), T);
   };

  cdouble VInv=0.0;
  double Sigma=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { VInv  += FA[Mu]*Y[Mu][Nu]*FB[Nu];
      Sigma += Theta*FA[Mu]*imag(Eps[Mu][Nu])*FB[Nu];
    };
  VInv  *= -1.0/ (Omega*Omega);
  Sigma *= 2.0*real(Omega) / (M_PI*ZVAC);
 
  I[0] = real(VInv);
  I[1] = imag(VInv);
  I[2] = Sigma; 
  
}

/***************************************************************/
/* For a given SWG basis function f_a, this routine computes   */
/* matrix elements of the VInverse and Sigma operators between */
/* f_a and all basis functions f_b that have nonzero overlap   */
/* with f_a (there are a maximum of MAXOVERLAP=7 such BFs.)    */
/* The indices of the f_b functions are returned in Indices.   */
/*                                                             */
/* VInverse = (1 - Eps^{-1}) / k^2                             */
/*                                                             */
/*    Sigma = Theta(T) * Im Eps * (2k / Pi*ZVAC)               */
/*                                                             */
/* If Temperature==NULL then Theta(T) is set to 1.             */
/***************************************************************/
int GetVInvAndSigmaEntries(SWGVolume *V,
                           int nfA, cdouble Omega,
                           IHAIMatProp *Temperature,
                           int Indices[MAXOVERLAP],
                           cdouble VInvEntries[MAXOVERLAP],
                           double SigmaEntries[MAXOVERLAP])
{
  Indices[0]=nfA;
  VInvEntries[0]=0.0;
  SigmaEntries[0]=0.0;
  int NNZ=1;

  SWGFace *FA = V->Faces[nfA];
  struct VIIData MyVIIData, *Data=&MyVIIData;
  Data->Omega       = Omega;
  Data->MP          = V->MP;
  Data->Temperature = Temperature;

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
        SigmaEntries[0] += I[2];
      }
     else
      { VInvEntries[NNZ]  = cdouble(I[0], I[1]);
        SigmaEntries[NNZ] = I[2];
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
        VInvEntries[0]  += cdouble( I[0], I[1] );
        SigmaEntries[0] += I[2];
      }
     else
      { VInvEntries[NNZ]  = cdouble( I[0], I[1] );
        SigmaEntries[NNZ] = I[2];
        Indices[NNZ]      = nfB;
        NNZ++;
      };

   }; // for(int ifB=0; ifB<4; ifB++)

  return NNZ;

} // routine GetVInvAndImEpsEntries

int GetVInvEntries(SWGVolume *V, int nfA, cdouble Omega,
                   int Indices[MAXOVERLAP],
                   cdouble VInvEntries[MAXOVERLAP])
{
  double SigmaEntries[MAXOVERLAP];
  return GetVInvAndSigmaEntries(V, nfA, Omega, 0,
                                Indices, VInvEntries, SigmaEntries);
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIBBICache *GCache = (noa==nob) ?  ObjectGCaches[noa] : 0;
  int InitialCacheSize = 0;
  if (GCache && GCache->Size()==0)
   { char CacheFileName[MAXSTR];
     snprintf(CacheFileName, MAXSTR, "%s.GCache",
              RemoveExtension(OA->MeshFileName));
     GCache->PreLoad(CacheFileName);
     InitialCacheSize = GCache->Size();
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Log("Assembling G(%i,%i)",noa,nob);
#ifndef USE_OPENMP
  Log(" no multithreading...");
#else
  int NumThreads=GetNumThreads();
  Log(" OpenMP multithreading (%i threads...)",NumThreads);
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
  if (GCache && GCache->Size()>InitialCacheSize)
   { char CacheFileName[MAXSTR];
     snprintf(CacheFileName, MAXSTR, "%s.GCache",
              RemoveExtension(OA->MeshFileName));
     GCache->Store(CacheFileName);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::AssembleVInvBlock(int no, cdouble Omega,
                                    IHAIMatProp *Temperature,
                                    SMatrix *VInv, SMatrix *Sigma,
                                    HMatrix *TInv,
                                    int Offset)
{
   Log("Adding VInv(%i)",no);
   SWGVolume *O = Objects[no];

#ifdef USE_OPENMP
   int NumThreads=GetNumThreads();
   if (VInv || Sigma)
    NumThreads=1;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
   for(int nr=0; nr<O->NumInteriorFaces; nr++)
    { int ncList[MAXOVERLAP];
      cdouble VInvEntries[MAXOVERLAP];
      double SigmaEntries[MAXOVERLAP];
      int NNZ = GetVInvAndSigmaEntries(O, nr, Omega, Temperature,
                                       ncList, VInvEntries, SigmaEntries);
      for(int nnz=0; nnz<NNZ; nnz++)
        { int nc = ncList[nnz];
          if (VInv)
           VInv->SetEntry(nr, nc, VInvEntries[nnz]);
          if (Sigma) 
           Sigma->SetEntry(nr, nc, SigmaEntries[nnz]);
          if (TInv)
           TInv->AddEntry(Offset+nr, Offset+nc, VInvEntries[nnz]);
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
      AssembleGBlock(noa, nob, Omega, M, 
                     BFIndexOffset[noa], BFIndexOffset[nob]);

      if (nob==noa)
       AssembleVInvBlock(noa, Omega, 0, 0, 0, M, BFIndexOffset[noa]);
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
