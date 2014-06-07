/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * TGMatrices.cc -- libSWG routines for forming the T and G matrix
 *                  blocks
 *
 * homer reid    -- 5/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libSGJC.h"
#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;

namespace buff {

#define MAXSTR 1000

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
#if 0
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
  
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetVInverseElement(SWGVolume *V, cdouble Omega, HMatrix *TInv)
{
#if 0
  for(int nfA=0; nfA<V->NumInteriorFaces; nfA++)
   { 
     SWGFace *FA = V->Faces[nfA];
     struct VIIData MyVIIData, *Data=&MyVIIData;
     Data->Omega = Omega;
     Data->MP    = V->MP;

     /*--------------------------------------------------------------*/
     /* handle interactions between face #nfA and face #nfB, where   */
     /* nfB runs over the 4 faces of the positive tetrahedron of nfA */
     /*--------------------------------------------------------------*/
     int ntA       = FA->iPTet;
     SWGTet *TA    = V->Tets[ntA];
     Data->QA      = V->Vertices + 3*(FA->iQP);
     Data->PreFacA = FA->Area / (3.0*TA->Volume);
     for(int ifB=0; ifB<4; ifB++)
      { 
        int nfB = TA->FI[ifB];
        if (nfB >= V->NumInteriorFaces) continue;

        SWGFace *FB = V->Faces[nfB];
        if ( FB->iPTet == ntA )
         { Data->QB      = V->Vertices + 3*(FB->iQP);
           Data->PreFacB = FB->Area / (3.0*TA->Volume);
         }
        else
         { Data->QB      = V->Vertices + 3*(FB->iQM);
           Data->PreFacB = -1.0*FB->Area / (3.0*TA->Volume);
         };

        cdouble Result, Error;
        TetInt(this, ntA, 0, 1.0, VInvIntegrand, (void *)Data,
               2, (double *)&Result, (double *)&Error, 0, 1.0e-4);
        M->AddEntry(nfA, nfB, Result);
      };

     /*--------------------------------------------------------------*/
     /* handle interactions between face #nfA and face #nfB, where   */
     /* nfB runs over the 3 faces of the negative tetrahedron of nfA */
     /* that differ from nfA                                         */
     /*--------------------------------------------------------------*/
     ntA           = FA->iMTet;
     TA            = V->Tets[ntA];
     Data->QA      = V->Vertices + 3*(FA->iQM);
     Data->PreFacA = -1.0*FA->Area / (3.0*TA->Volume);
     for(int ifB=0; ifB<4; ifB++)
      { 
        int nfB = TA->FI[ifB];
        if (nfB >= V-NumInteriorFaces) continue;
        if (nfB==nfA) continue;

        SWGFace *FB = V->Faces[nfB];
        if ( FB->iPTet == ntA )
         { Data->QB      = V->Vertices + 3*(FB->iQP);
           Data->PreFacB = FB->Area / (3.0*TA->Volume);
         }
        else
         { Data->QB      = V->Vertices + 3*(FB->iQM);
           Data->PreFacB = -1.0*FB->Area / (3.0*TA->Volume);
         };

        cdouble Result, Error;
        TetInt(this, ntA, 0, 1.0, VInvIntegrandIntegrand, (void *)Data,
               2, (double *)&Result, (double *)&Error, 0, 1.0e-4);
        M->AddEntry(nfA, nfB, Result);
      };

   }; // for(int nfA=0; nfA<NumInteriorFaces; nfA++)

#endif
} // routine AddPotentialTermToTInverse

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GSurfaceIntegrand(double *xA, double *bA,
                       double DivbA, double *nHatA,
                       double *xB, double *bB,
                       double DivbB, double *nHatB,
                       void *UserData, double *I)
{
#if 0
  double R[3];
  R[0] = (xA[0]-xB[0]);
  R[1] = (xA[1]-xB[1]);
  R[2] = (xA[2]-xB[2]);
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r = sqrt(r2);

  IntegrandData *Data = (IntegrandData *)UserData;
  cdouble k = Data->k;

  cdouble Xi = II*k*r;
  cdouble h, p, w;
  Gethpq(Xi, &h, &p, &q);

  double NdotN =   nHatA[0]*nHatB[0] 
                 + nHatA[1]*nHatB[1] 
                 + nHatA[2]*nHatB[2];

  double PhiDotPhi = bA[0]*bB[0] + bA[1]*bB[1] * bA[2]*bB[2];

  double VmVDotXi  = bA[0]*bB[0] + bA[1]*bB[1] * bA[2]*bB[2];
  
  cdouble *zI = (cdouble *)I;
  //zI[0] = -1.0 * DotProd * h / (II*k);

#endif
}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions using the surface-integral approach of Bleszynski */
/* et. al.                                                     */
/***************************************************************/
cdouble GetUMatrixElement_SI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega)
{
#if 0
  cdouble RetVal=0.0;
  int fDim=2;

  // 36 surface integrals
  for(int ASign=+1; ASign>=-1; ASign-=2)
   for(int BSign=+1; BSign>=-1; BSign-=2)
    {
      int ntA    = (ASign==1) ? FA->iPTet : FA->iMTet;
      int iQA    = (ASign==1) ? FA->iQP   : FA->iQM;
      Data->QA   = VA->Vertices + 3*iQA;
      SWGTet *TA = VA->Tets[ ntA ];

      int ntB = (BSign==1) ? FB->iPTet : FB->iMTet;
      int iQB = (BSign==1) ? FB->iQP   : FB->iQM;
      SWGTet *TB = VB->Tets[ ntB ];

      for(int nfP=0; nfP<4; nfP++)
       for(int nfQ=0; nfQ<4; nfQ++)
        { 
          if ( (TA->FI[nfP] == nfA) || (TB->FI[nfQ] == nfB) )
           continue;

          double PResult[2], PError[2];
          FaceFaceInt(V, ntA, nfP, iQA, ASign,
                      V, ntB, nfQ, iQB, BSign,
                      GSurfaceIntegrand, (void *)Data, fDim,
                      PResult, PError, 100000, 1.0e-8); 
          RetVal += cdouble( PResult[0], PResult[1] );
        };
    };

  return RetVal;

#endif
}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions in the dipole approximation retaining NumTerms    */
/* terms (here NumTerms may be 1 or 2).                        */
/***************************************************************/
cdouble GetUMatrixElement_DA(SWGVolume *VA, int nfA, 
                             SWGVolume *VB, int nfB,
                             cdouble Omega, int NumTerms)
{}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeGMatrix(SWGVolume *VA, SWGVolume *VB,
                    cdouble Omega, HMatrix *G)
{
#if 0
  for(int nfA=0; nfA<VA->NumInteriorFaces; nfA++)
   for(int nfB=0; nfB<VB->NumInteriorFaces; nfB++)
    G->SetEntry(nfA, nfB, ComputeGMatrixEntry(VA, nfA, VB, nfB, Omega));
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AssembleVIEMatrixBlock(int noa, int nob, cdouble Omega,
                                             HMatrix *M, int RowOffset, int ColOffset)
{
#if 0
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGVolume *OA = Objects[noa];
  SWGVolume *OB = Objects[nob];
  int SameObject = (noa==nob) ? 1 : 0;

  for(int nfa=0; nfa<OA->NumInteriorFaces; nfa++)
   for(int nfb=SameObject*nfa; nfb<OB->NumInteriorFaces; nfb++)
    M->SetEntry(RowOffset + nfa, ColOffset + nfb,
                GetGMatrixElement(OA, nfa, OB, nfb, Omega) );

#endif
  /***************************************************************/
  /**************************************************************/
  /***************************************************************/
#if 0
  if (noa==nob)
   { 
     for(int nfa=0; nfa<OA->NumInteriorFaces; nfa++)
      GetTInverseMatrixEntries(
   };
#endif

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AssembleVIEMatrix(cdouble Omega, HMatrix *TInv)
{
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::AllocateVIEMatrix()
{
  return new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *SWGGeometry::AllocateRHSVector()
{
  return new HVector(TotalBFs, LHM_COMPLEX);
}

} // namespace buff
