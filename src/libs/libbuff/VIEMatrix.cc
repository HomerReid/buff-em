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

T1=0.0;
T2=-DivbA*DivbB*h/k2;

  zI[0] = NdotN * (T1 + T2) / (-4.0*M_PI*II*k);

  //zI[0] = NdotN * h / (-4.0*M_PI*II*k);
//printf("T1, T2 = (%e,%e)\n", 
//        abs(NdotN * T1 / (-4.0*M_PI*II*k)),
//        abs(NdotN * T2 / (-4.0*M_PI*II*k)));
//printf("(DivbA*DivbB)/(9*k2): %s \n",CD2S(DivbA*DivbB / (9.0*k2)));
//printf("(h,w,p,k2*DQdotR*p)=(%e,%e,%e,%e)\n",
//        abs(h),abs(w),abs(p),abs(k2*DQdotR*p));

}

/***************************************************************/
/* Compute the G-matrix element between two tetrahedral basis  */
/* functions using the surface-integral approach of Bleszynski */
/* et. al.                                                     */
/***************************************************************/
cdouble GetGMatrixElement_SI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega)
{
  int fDim=2;

  GSIData MyData, *Data = &MyData;

  Data->k = Omega;
  
  SWGFace *FA = VA->Faces[nfA];
  SWGFace *FB = VB->Faces[nfB];

  // 36 surface integrals
  cdouble RetVal=0.0;
  for(int ASign=+1; ASign>=-1; ASign-=2)
   for(int BSign=+1; BSign>=-1; BSign-=2)
    {
      int ntA    = (ASign==1) ? FA->iPTet : FA->iMTet;
      int iQA    = (ASign==1) ? FA->iQP   : FA->iQM;
      Data->QA   = VA->Vertices + 3*iQA;
      SWGTet *TA = VA->Tets[ ntA ];

      int ntB    = (BSign==1) ? FB->iPTet : FB->iMTet;
      int iQB    = (BSign==1) ? FB->iQP   : FB->iQM;
      Data->QB   = VB->Vertices + 3*iQB;
      SWGTet *TB = VB->Tets[ ntB ];

      for(int nfP=0; nfP<4; nfP++)
       for(int nfQ=0; nfQ<4; nfQ++)
        { 
          if ( (TA->FI[nfP] == nfA) || (TB->FI[nfQ] == nfB) )
           continue;

          double PResult[2], PError[2];
          FaceFaceInt(VA, ntA, nfP, iQA, ASign,
                      VB, ntB, nfQ, iQB, BSign,
                      GSurfaceIntegrand, (void *)Data, fDim,
                      PResult, PError, 0, 1000, 1.0e-8);
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
                  int NeedQ=true)
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
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#if 0
double R[3];
R[0] = FA->Centroid[0] - FB->Centroid[0];
R[1] = FA->Centroid[1] - FB->Centroid[1];
R[2] = FA->Centroid[2] - FB->Centroid[2];
double r=sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );
return (OA->Tets[FA->iPTet]->Volume + OA->Tets[FA->iMTet]->Volume)
      *(OB->Tets[FB->iPTet]->Volume + OB->Tets[FB->iMTet]->Volume)
      *exp(II*Omega*r) / (4.0*M_PI*r);
#endif
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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

  double *x0A = FA->Centroid;
  double *x0B = FB->Centroid;

  double R[3];
  R[0] = x0A[0] - x0B[0];
  R[1] = x0A[1] - x0B[1];
  R[2] = x0A[2] - x0B[2];
  double r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);

  if ( r > 20.0*fmax( FA->Radius, FB->Radius) )
   return GetGMatrixElement_DA(VA, nfA, VB, nfB, Omega, 1);
  else if ( r > 10.0*fmax( FA->Radius, FB->Radius) )
   return GetGMatrixElement_DA(VA, nfA, VB, nfB, Omega, 2);
  else 
   return GetGMatrixElement_SI(VA, nfA, VB, nfB, Omega);
  
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
  int SameObject = (noa==nob) ? 1 : 0;

  for(int nfa=0; nfa<OA->NumInteriorFaces; nfa++)
   for(int nfb=SameObject*nfa; nfb<OB->NumInteriorFaces; nfb++)
    M->SetEntry(RowOffset + nfa, ColOffset + nfb,
                GetGMatrixElement(OA, nfa, OB, nfb, Omega) );

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
