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
 * Cubature.cc -- libSWG routines for integrals over tetrahedra
 *                and triangular faces
 *
 * homer reid  -- 5/2014
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
typedef void (*UserTIntegrand)(double *x, double *b, double Divb,
                               void *UserData, double *I);

typedef void (*UserTTIntegrand)(double *xA, double *bA, double DivbA,
                                double *xB, double *bB, double DivbB,
                                void *UserData, double *I);

typedef void (*UserFIntegrand)(double *x, double *b, double Divb, double *nHat,
                               void *UserData, double *I);

typedef void (*UserFFIntegrand)(double *xA, double *bA, double DivbA, double *nHatA,
                                double *xB, double *bB, double DivbB, double *nHatB,
                                void *UserData, double *I);

/***************************************************************/
/* utility routine to get the outward-pointing normal to a     */
/* triangle with vertices (V1, V2, V3). 'Outward-pointing'     */
/* means 'pointing away from RefPnt.'                          */
/***************************************************************/
void GetOutwardPointingNormal(double *V1, double *V2, double *V3, 
                              double *RefPnt, double *nHat)
{
  double L1[3], L2[3]; 
  VecSub(V2, V1, L1);
  VecSub(V3, V1, L2);
  VecCross(L1, L2, nHat);
  VecNormalize(nHat);

  double Centroid[3];
  Centroid[0] = (V1[0] + V2[0] + V3[0])/3.0;
  Centroid[1] = (V1[1] + V2[1] + V3[1])/3.0;
  Centroid[2] = (V1[2] + V2[2] + V3[2])/3.0;

  double D1 = VecDistance(Centroid, RefPnt);

  VecPlusEquals(Centroid, 0.1*VecNorm(L1), nHat );
  double D2 = VecDistance(Centroid, RefPnt);

  if (D2<D1)
   VecScale(nHat, -1.0);
}

/*=============================================================*/
/*= PART 1: integrals over single tetrahedra ==================*/
/*=============================================================*/
/***************************************************************/
/* internal data structure passed to TIIntegrand ***************/
/***************************************************************/  
typedef struct TIData 
 { 
   double *Q, L1[3], L2[3], L3[3]; 
   double Volume, PreFac;
   void *UserData;
   UserTIntegrand Integrand;

 } TIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int TIIntegrand(unsigned ndim, const double *uvw, void *params,
                unsigned fdim, double *fval)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TIData *Data   = (TIData *)params;
  double *Q      = Data->Q;
  double *L1     = Data->L1;
  double *L2     = Data->L2;
  double *L3     = Data->L3;
  double Volume  = Data->Volume;
  double PreFac  = Data->PreFac;
  void *UserData = Data->UserData;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double u  = uvw[0];
  double v  = uvw[1];
  double w  = uvw[2];
  double vp = (1.0-u)*v;
  double wp = (1.0-u)*(1.0-v)*w;
  double Jacobian = (1.0-u)*(1.0-u)*(1.0-v) * 6.0 * Volume;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double x[3], b[3];
  for(int Mu=0; Mu<3; Mu++)
   { 
     b[Mu] = u*L1[Mu] + vp*L2[Mu] + wp*L3[Mu];
     x[Mu] = Q[Mu] + b[Mu];
     b[Mu] *= PreFac;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Data->Integrand(x, b, 3.0*PreFac, UserData, fval);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int n=0; n<fdim; n++)
   fval[n]*=Jacobian;

  return 0;
  
}


/***************************************************************/
/* evaluate an integral over the volume of a single tetrahedron*/
/***************************************************************/
void TetInt(SWGVolume *V, int nt, int iQ, double Sign,
             UserTIntegrand Integrand, void *UserData,
            int fdim, double *Result, double *Error,
            int MaxEvals, double RelTol)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGTet *T     = V->Tets[nt];
  double *Q     = V->Vertices + 3*(T->VI[ iQ ]);
  double *V1    = V->Vertices + 3*(T->VI[ (iQ+1)%4 ]);
  double *V2    = V->Vertices + 3*(T->VI[ (iQ+2)%4 ]);
  double *V3    = V->Vertices + 3*(T->VI[ (iQ+3)%4 ]);
  double PreFac = Sign*V->Faces[T->FI[iQ]]->Area / (3.0 * T->Volume);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TIData MyTIData, *Data= &(MyTIData);
  Data->Q = Q;
  VecSub(V1, Q, Data->L1);
  VecSub(V2, Q, Data->L2);
  VecSub(V3, Q, Data->L3);
  Data->Volume    = T->Volume;
  Data->PreFac    = PreFac;
  Data->UserData  = UserData;
  Data->Integrand = Integrand;

  double Lower[3]={0.0, 0.0, 0.0};
  double Upper[3]={1.0, 1.0, 1.0};
  hcubature(fdim, TIIntegrand, (void *)Data, 3, Lower, Upper,
	    MaxEvals, 0.0, RelTol, ERROR_INDIVIDUAL, Result, Error);

}

/*=============================================================*/
/*= PART 2: integrals over pairs of tetrahedra ================*/
/*=============================================================*/

/***************************************************************/
/* internal data structure passed to TTIIntegrand **************/
/***************************************************************/  
typedef struct TTIData 
 { 
   double *QA, L1A[3], L2A[3], L3A[3]; 
   double VolumeA, PreFacA;

   double *QB, L1B[3], L2B[3], L3B[3]; 
   double VolumeB, PreFacB;

   void *UserData;
   UserTTIntegrand Integrand;

 } TTIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int TTIIntegrand(unsigned ndim, const double *uvw, void *params,
                 unsigned fdim, double *fval)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TTIData *Data= (TTIData *)params;

  double *QA     = Data->QA;
  double *L1A    = Data->L1A;
  double *L2A    = Data->L2A;
  double *L3A    = Data->L3A;
  double VolumeA = Data->VolumeA;
  double PreFacA = Data->PreFacA;

  double *QB     = Data->QB;
  double *L1B    = Data->L1B;
  double *L2B    = Data->L2B;
  double *L3B    = Data->L3B;
  double VolumeB = Data->VolumeB;
  double PreFacB = Data->PreFacB;

  void *UserData = Data->UserData;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double uA = uvw[0];
  double vA = uvw[1];
  double wA = uvw[2];
  double vpA = (1.0-uA)*vA;
  double wpA = (1.0-uA)*(1.0-vA)*wA;
  double JacobianA = (1.0-uA)*(1.0-uA)*(1.0-vA) * 6.0 * VolumeA;

  double uB = uvw[3];
  double vB = uvw[4];
  double wB = uvw[5];
  double vpB = (1.0-uB)*vB;
  double wpB = (1.0-uB)*(1.0-vB)*wB;
  double JacobianB = (1.0-uB)*(1.0-uB)*(1.0-vB) * 6.0 * VolumeB;

  double Jacobian = JacobianA * JacobianB;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double xA[3], bA[3], xB[3], bB[3];
  for(int Mu=0; Mu<3; Mu++)
   { 
     bA[Mu] = uA*L1A[Mu] + vpA*L2A[Mu] + wpA*L3A[Mu];
     xA[Mu] = QA[Mu] + bA[Mu];
     bA[Mu] *= PreFacA;

     bB[Mu] = uB*L1B[Mu] + vpB*L2B[Mu] + wpB*L3B[Mu];
     xB[Mu] = QB[Mu] + bB[Mu];
     bB[Mu] *= PreFacB;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Data->Integrand(xA, bA, 3.0*PreFacA, xB, bB, 3.0*PreFacB, UserData, fval);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int n=0; n<fdim; n++)
   fval[n]*=Jacobian;

  return 0;
  
}


/***************************************************************/
/* evaluate an integral over a pair of tetrahedra              */
/***************************************************************/
void TetTetInt(SWGVolume *VA, int ntA, int iQA, double SignA,
               SWGVolume *VB, int ntB, int iQB, double SignB,
               UserTTIntegrand Integrand, void *UserData,
               int fdim, double *Result, double *Error,
               int MaxEvals, double RelTol)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGTet *TA     = VA->Tets[ntA];
  double *QA     = VA->Vertices + 3*(TA->VI[ iQA ]);
  double *V1A    = VA->Vertices + 3*(TA->VI[ (iQA+1)%4 ]);
  double *V2A    = VA->Vertices + 3*(TA->VI[ (iQA+2)%4 ]);
  double *V3A    = VA->Vertices + 3*(TA->VI[ (iQA+3)%4 ]);
  double PreFacA = SignA*(VA->Faces[TA->FI[iQA]]->Area) / (3.0 * TA->Volume);

  SWGTet *TB     = VB->Tets[ntB];
  double *QB     = VB->Vertices + 3*(TB->VI[ iQB ]);
  double *V1B    = VB->Vertices + 3*(TB->VI[ (iQB+1)%4 ]);
  double *V2B    = VB->Vertices + 3*(TB->VI[ (iQB+2)%4 ]);
  double *V3B    = VB->Vertices + 3*(TB->VI[ (iQB+3)%4 ]);
  double PreFacB = SignB*(VB->Faces[TB->FI[iQB]]->Area) / (3.0 * TB->Volume);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TTIData MyTTIData, *Data= &(MyTTIData);

  Data->QA = QA;
  VecSub(V1A, QA, Data->L1A);
  VecSub(V2A, QA, Data->L2A);
  VecSub(V3A, QA, Data->L3A);
  Data->VolumeA  = TA->Volume;
  Data->PreFacA  = PreFacA;

  Data->QB = QB;
  VecSub(V1B, QB, Data->L1B);
  VecSub(V2B, QB, Data->L2B);
  VecSub(V3B, QB, Data->L3B);
  Data->VolumeB  = TB->Volume;
  Data->PreFacB  = PreFacB;

  Data->UserData  = UserData;
  Data->Integrand = Integrand;

  double Lower[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Upper[6]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  hcubature(fdim, TTIIntegrand, (void *)Data, 6, Lower, Upper,
	    MaxEvals, 0.0, RelTol, ERROR_INDIVIDUAL, Result, Error);

}

/*=============================================================*/
/*= PART 3: integrals over single faces =======================*/
/*=============================================================*/

/***************************************************************/
/* internal data structure passed to FIIntegrand ***************/
/***************************************************************/  
typedef struct FIData 
 { 
   double *V1, L1[3], L2[3], *Q, nHat[3];
   double Area, PreFac;
   void *UserData;
   UserFIntegrand Integrand;

 } FIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int FIIntegrand(unsigned ndim, const double *uv, void *params,
                unsigned fdim, double *fval)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIData *Data   = (FIData *)params;
  double *V1     = Data->V1;
  double *L1     = Data->L1;
  double *L2     = Data->L2;
  double *Q      = Data->Q;
  double *nHat   = Data->nHat;
  double Area    = Data->Area;  
  double PreFac  = Data->PreFac;
  void *UserData = Data->UserData;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double u  = uv[0];
  double v  = uv[1];
  double vp = (1.0-u)*v;
  double Jacobian = (1.0-u) * 2.0 * Area;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double x[3], b[3];
  for(int Mu=0; Mu<3; Mu++)
   { 
     x[Mu] = V1[Mu] + u*L1[Mu] + vp*L2[Mu];
     b[Mu] = PreFac * (x[Mu] - Q[Mu]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Data->Integrand(x, b, 3.0*PreFac, nHat, UserData, fval);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int n=0; n<fdim; n++)
   fval[n]*=Jacobian;

  return 0;
  
}


/***************************************************************/
/* evaluate an integral over the area of a single face         */
/***************************************************************/
void FaceInt(SWGVolume *V, int nt, int nf, int iQ, double Sign,
             UserFIntegrand Integrand, void *UserData,
             int fdim, double *Result, double *Error,
             int MaxEvals, double RelTol)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGTet  *T    = V->Tets[nt];
  SWGFace *F    = V->Faces[ T->FI[nf] ];
  double *V1    = V->Vertices + 3*(F->iV1);
  double *V2    = V->Vertices + 3*(F->iV2);
  double *V3    = V->Vertices + 3*(F->iV3);
  double *Q     = V->Vertices + 3*iQ;
  double PreFac = Sign*F->Area / (3.0*T->Volume);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIData MyFIData, *Data= &(MyFIData);
  VecSub(V2, V1, Data->L1);
  VecSub(V3, V1, Data->L2);
  Data->V1        = V1;
  Data->Q         = Q;
  Data->Area      = F->Area;
  Data->PreFac    = PreFac;
  Data->UserData  = UserData;
  Data->Integrand = Integrand;
  GetOutwardPointingNormal(V1, V2, V3, T->Centroid, Data->nHat);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  hcubature(fdim, FIIntegrand, (void *)Data, 2, Lower, Upper,
	    MaxEvals, 0.0, RelTol, ERROR_INDIVIDUAL, Result, Error);

}

/*=============================================================*/
/*= PART 4: integrals over pairs of faces. ====================*/
/*=============================================================*/

/***************************************************************/
/* internal data structure passed to FFIIntegrand **************/
/***************************************************************/  
typedef struct FFIData 
 { 
   double *V1A, L1A[3], L2A[3], *QA, nHatA[3];
   double AreaA, PreFacA;

   double *V1B, L1B[3], L2B[3], *QB, nHatB[3];
   double AreaB, PreFacB;

   void *UserData;
   UserFFIntegrand Integrand;

 } FFIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int FFIIntegrand(unsigned ndim, const double *uv, void *params,
                 unsigned fdim, double *fval)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FFIData *Data   = (FFIData *)params;

  double *V1A     = Data->V1A;
  double *L1A     = Data->L1A;
  double *L2A     = Data->L2A;
  double *QA      = Data->QA;
  double *nHatA   = Data->nHatA;
  double AreaA    = Data->AreaA;
  double PreFacA  = Data->PreFacA;

  double *V1B     = Data->V1B;
  double *L1B     = Data->L1B;
  double *L2B     = Data->L2B;
  double *QB      = Data->QB;
  double *nHatB   = Data->nHatB;
  double AreaB    = Data->AreaB;
  double PreFacB  = Data->PreFacB;

  void *UserData = Data->UserData;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double uA  = uv[0];
  double vA  = uv[1];
  double vpA = (1.0-uA)*vA;
  double JacobianA = (1.0-uA) * 2.0 * AreaA;

  double uB  = uv[2];
  double vB  = uv[3];
  double vpB = (1.0-uB)*vB;
  double JacobianB = (1.0-uB) * 2.0 * AreaB;

  double Jacobian = JacobianA * JacobianB;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double xA[3], bA[3], xB[3], bB[3];
  for(int Mu=0; Mu<3; Mu++)
   { 
     xA[Mu] = V1A[Mu] + uA*L1A[Mu] + vpA*L2A[Mu];
     bA[Mu] = PreFacA * (xA[Mu] - QA[Mu]);

     xB[Mu] = V1B[Mu] + uB*L1B[Mu] + vpB*L2B[Mu];
     bB[Mu] = PreFacB * (xB[Mu] - QB[Mu]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Data->Integrand(xA, bA, 3.0*PreFacA, nHatA,
                  xB, bB, 3.0*PreFacB, nHatB, UserData, fval);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int n=0; n<fdim; n++)
   fval[n]*=Jacobian;

  return 0;
  
}

/***************************************************************/
/* evaluate an integral over two faces                         */
/***************************************************************/
void FaceFaceInt(SWGVolume *VA, int ntA, int nfA, int iQA, double SignA,
                 SWGVolume *VB, int ntB, int nfB, int iQB, double SignB,
                 UserFFIntegrand Integrand, void *UserData,
                 int fdim, double *Result, double *Error,
                 int MaxEvals, double RelTol)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGTet  *TA    = VA->Tets[ntA];
  SWGFace *FA    = VA->Faces[ TA->FI[nfA] ];
  double *V1A    = VA->Vertices + 3*(FA->iV1);
  double *V2A    = VA->Vertices + 3*(FA->iV2);
  double *V3A    = VA->Vertices + 3*(FA->iV3);
  double  *QA    = VA->Vertices + 3*iQA;
  double PreFacA = SignA * FA->Area / (3.0*TA->Volume);

  SWGTet  *TB    = VB->Tets[ntB];
  SWGFace *FB    = VB->Faces[ TB->FI[nfB] ];
  double *V1B    = VB->Vertices + 3*(FB->iV1);
  double *V2B    = VB->Vertices + 3*(FB->iV2);
  double *V3B    = VB->Vertices + 3*(FB->iV3);
  double  *QB    = VB->Vertices + 3*iQB;
  double PreFacB = SignB * FB->Area / (3.0*TB->Volume);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FFIData MyFIData, *Data= &(MyFIData);

  Data->V1A = V1A;
  VecSub(V2A, V1A, Data->L1A);
  VecSub(V3A, V1A, Data->L2A);
  Data->QA  = QA;
  Data->AreaA     = FA->Area;
  Data->PreFacA   = PreFacA;
  GetOutwardPointingNormal(V1A, V2A, V3A, TA->Centroid, Data->nHatA);

  Data->V1B = V1B;
  VecSub(V2B, V1B, Data->L1B);
  VecSub(V3B, V1B, Data->L2B);
  Data->QB  = QB;
  Data->AreaB     = FB->Area;
  Data->PreFacB   = PreFacB;
  GetOutwardPointingNormal(V1B, V2B, V3B, TB->Centroid, Data->nHatB);

  Data->UserData  = UserData;
  Data->Integrand = Integrand;

  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  hcubature(fdim, FFIIntegrand, (void *)Data, 4, Lower, Upper,
	    MaxEvals, 0.0, RelTol, ERROR_INDIVIDUAL, Result, Error);

}

} // namespace buff
