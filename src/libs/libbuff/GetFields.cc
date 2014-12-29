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
 * GetFields.cc  -- libbuff class methods for computing scattered
 *               -- electric and magnetic fields 
 *
 * homer reid    -- 3/2007  -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libMDInterp.h>

#include "libscuff.h"
#include "libbuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define MAXFUNC 50
#define II cdouble(0,1)

using namespace scuff;

namespace scuff{

void CalcGC(double R1[3], double R2[3],
            cdouble Omega, cdouble EpsR, cdouble MuR, 
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

               }

namespace buff {

/***************************************************************/
/* get 1BF fields using surface-integral method ****************/
/***************************************************************/
#if 0

void ExpRel23(cdouble x, cdouble *ExpRel2, cdouble *ExpRel3);

typedef struct G1BFData
 {
   cdouble k;
   double *XEval;
 } G1BFData;

void G1BFSIIntegrand(double *XSource, double *b, double Divb, double *nHat,
                     void *UserData, double *I)
{ 
  /*--------------------------------------------------------------*
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  G1BFData *Data = (G1BFData *)UserData;
  cdouble k     = Data->k;
  double *XEval = Data->XEval;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double R[3];
  R[0] = XEval[0] - XSource[0];
  R[1] = XEval[1] - XSource[1];
  R[2] = XEval[2] - XSource[2];
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r = sqrt(r2);

  /*--------------------------------------------------------------*/
  /*- i think the following calculation could probably be         */
  /*- accelerated -- in particular, do I really need to be making */
  /*- two separate calls to ExpRelBar()?                          */
  /*--------------------------------------------------------------*/
  cdouble Xi  = II*k*r; 
  cdouble Xi2 = Xi*Xi;
  cdouble Xi4 = Xi2*Xi2;
  cdouble ER2, ER3;
  ExpRel23(Xi, &ER2, &ER3);
  cdouble h = ER2 / (4.0*M_PI*Xi);
  ExpRel23(-Xi, &ER2, &ER3);
  cdouble ExpFac = exp(II*k*r) / (4.0*M_PI*Xi);
  cdouble q = ExpFac * ER2 / Xi2;
  cdouble t = -1.0 * ExpFac * (ER2 + 2.0*ER3) / Xi4;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double rdnHat = R[0]*nHat[0] + R[1]*nHat[1] + R[2]*nHat[2];

  cdouble k2     = k*k;
  double PreFac  = ZVAC*Divb/3.0;
  cdouble Coeff1 = ZVAC * k2 * rdnHat * q;
  cdouble Coeff2 = PreFac * ( 3.0*q - h ); 
  cdouble Coeff3 = PreFac * 3.0*k2*rdnHat*t;

  cdouble *zI = (cdouble *)I;
  zI[0] = Coeff1*b[0] + Coeff2*nHat[0] - Coeff3*R[0];
  zI[1] = Coeff1*b[1] + Coeff2*nHat[1] - Coeff3*R[1];
  zI[2] = Coeff1*b[2] + Coeff2*nHat[2] - Coeff3*R[2];
  zI[3] = 0.0;
  zI[4] = 0.0;
  zI[5] = 0.0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Get1BFFields_SI(SWGVolume *O, int nf, cdouble k, double X[3],
                     cdouble EH[6], int Order=0)
{
  SWGFace *F=O->Faces[nf];

  G1BFData MyData, *Data=&MyData;
  Data->k = k;
  Data->XEval = X;
  
  EH[0]=EH[1]=EH[2]=EH[3]=EH[4]=EH[5]=0.0;
  cdouble pResult[6], pError[6];
  for(int Sign=1; Sign>=-1; Sign-=2)
   for(int iF=0; iF<4; iF++)
    { 
      int nT   = (Sign==1) ? F->iPTet  : F->iMTet;
      int iFBF = (Sign==1) ? F->PIndex : F->MIndex;

      FaceInt(O, nT, iF, iFBF, Sign, G1BFSIIntegrand, (void *)Data,
              12, (double *)pResult, (double *)pError, Order, 0, 1.0e-4);
      EH[0] += pResult[0];
      EH[1] += pResult[1];
      EH[2] += pResult[2];
      EH[3] += pResult[3];
      EH[4] += pResult[4];
      EH[5] += pResult[5];
    };

}
#endif

/***************************************************************/
/* get 1BF fields using volume-integral method *****************/
/***************************************************************/
typedef struct GFVIData
 {
   cdouble Omega;
   double *XDest;

 } GFVIData;

void GetFields_VIntegrand(double *XSource, double *b, double Divb,
                          void *UserData, double *I)
{
  GFVIData *Data = (GFVIData *)UserData;
  cdouble Omega = Data->Omega;
  double *XDest = Data->XDest;
   
  cdouble GMuNu[3][3], CMuNu[3][3];
  CalcGC(XDest, XSource, Omega, 1.0, 1.0, GMuNu, CMuNu, 0, 0);

  cdouble EPreFac = II*Omega*ZVAC;
  cdouble HPreFac = -II*Omega;
  cdouble *zI =  (cdouble *)I;
  zI[0] = EPreFac*(GMuNu[0][0]*b[0] + GMuNu[0][1]*b[1] + GMuNu[0][2]*b[2]);
  zI[1] = EPreFac*(GMuNu[1][0]*b[0] + GMuNu[1][1]*b[1] + GMuNu[1][2]*b[2]);
  zI[2] = EPreFac*(GMuNu[2][0]*b[0] + GMuNu[2][1]*b[1] + GMuNu[2][2]*b[2]);
  zI[3] = HPreFac*(CMuNu[0][0]*b[0] + CMuNu[0][1]*b[1] + CMuNu[0][2]*b[2]);
  zI[4] = HPreFac*(CMuNu[1][0]*b[0] + CMuNu[1][1]*b[1] + CMuNu[1][2]*b[2]);
  zI[5] = HPreFac*(CMuNu[2][0]*b[0] + CMuNu[2][1]*b[1] + CMuNu[2][2]*b[2]);

  CalcGC(XDest, XSource, Omega, 1.0, 1.0, GMuNu, CMuNu, 0, 0);
  
}

void Get1BFFields_VI(SWGVolume *O, int nf, cdouble Omega, double X[3],
                     cdouble EH[6], int NumPts=0)
{
  GFVIData MyData, *Data=&MyData;
  Data->Omega = Omega;
  Data->XDest = X;

  double Error[12];
  BFInt(O, nf, GetFields_VIntegrand, (void *)Data,
        12, (double *)EH, Error, NumPts, 0, 1.0e-6);

}

/***************************************************************/
/* get 1BF fields using dipole / quadrupole approximation      */
/***************************************************************/
#if 0
void Get1BFFields_DA(SWGVolume *O, int nf, cdouble Omega, double X[3],
                     cdouble EH[6], int NumTerms=1)
{
  SWGFace *F = O->Faces[nf];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double J[3], Q[3][3];
  cdouble GMuNu[3][3], GMuNuRho[3][3][3];
  cdouble CMuNu[3][3], CMuNuRho[3][3][3];

  GetDQMoments(O, nf, J, Q, (NumTerms==1) ? false : true );
  CalcGC(X, F->Centroid, Omega, 1.0, 1.0, GMuNu, CMuNu,
         (NumTerms==1) ? 0 : GMuNuRho, (NumTerms==1) ? 0 : CMuNuRho);
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  EH[0]=EH[1]=EH[2]=EH[3]=EH[4]=EH[5]=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { EH[ Mu + 0 ] += GMuNu[Mu][Nu]*J[Nu];
      EH[ Mu + 3 ] += CMuNu[Mu][Nu]*J[Nu];
    };

  if (NumTerms>1)
   { 
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       for(int Rho=0; Rho<3; Rho++)
        { EH[ Mu + 0 ] -= GMuNuRho[Mu][Nu][Rho]*Q[Nu][Rho];
          EH[ Mu + 3 ] -= CMuNuRho[Mu][Nu][Rho]*Q[Nu][Rho];
        };
   };

  cdouble EPreFac = II*Omega*ZVAC, HPreFac = -II*Omega;
  EH[0] *= EPreFac;
  EH[1] *= EPreFac;
  EH[2] *= EPreFac;
  EH[3] *= HPreFac;
  EH[4] *= HPreFac;
  EH[5] *= HPreFac;

}
#endif

/***************************************************************/
/* get the E and H fields due to a single SWG basis function   */
/* populated with unit strength                                */
/***************************************************************/
void Get1BFFields(SWGVolume *O, int nf, cdouble Omega, double X[3],
                  cdouble EH[6])
{
  SWGFace *F = O->Faces[nf];
  double rRel = VecDistance(X, F->Centroid) / F->Radius;
  if (rRel>10.0)
   Get1BFFields_VI(O, nf, Omega, X, EH, 4);
  else //if (rRel>1.0)
   Get1BFFields_VI(O, nf, Omega, X, EH, 16);
#if 0
  else
   Get1BFFields_SI(O, nf, Omega, X, EH, 20);
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::GetFields(IncField *IF, HVector *J, cdouble Omega,
                                HMatrix *XMatrix, HMatrix *FMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( XMatrix==0 || XMatrix->NC!=3 || XMatrix->NR==0 )
   ErrExit("wrong-size XMatrix (%ix%i) passed to GetFields",XMatrix->NR,XMatrix->NC);

  if (FMatrix==0) 
   FMatrix=new HMatrix(XMatrix->NR, 6, LHM_COMPLEX);
  else if ( (FMatrix->NR != XMatrix->NR) || (FMatrix->NC!=6) ) 
   { Warn(" ** warning: wrong-size FMatrix passed to GetFields(); allocating new matrix");
     FMatrix=new HMatrix(XMatrix->NR, 6, LHM_COMPLEX);
   };
  FMatrix->Zero();

  if (IF)
   IF->SetFrequency(Omega, true);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  // real and imaginary parts of the E and H field components,
  // which (as far as I can tell) need to be double-valued 
  // scalars to use with the OpenMP 'reduction' clause below...
  // if anybody knows how to use a cdouble-valued array in 
  // a 'reduction' clause please let me know!
  double ExReal, ExImag, EyReal, EyImag, EzReal, EzImag;
  double HxReal, HxImag, HyReal, HyImag, HzReal, HzImag;

  int NumTasks, NumThreads = GetNumThreads();
  Log("Getting fields...");
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:ExReal, ExImag, EyReal, EyImag, EzReal, EzImag, \
                                     HxReal, HxImag, HyReal, HyImag, HzReal, HzImag)
#endif
  for(int nr=0; nr<XMatrix->NR; nr++)
   { 
     double X[3];
     X[0] = XMatrix->GetEntryD(nr, 0);
     X[1] = XMatrix->GetEntryD(nr, 1);
     X[2] = XMatrix->GetEntryD(nr, 2);

     cdouble EHInc[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
     if (IF) IF->GetFields(X, EHInc);

     ExReal=ExImag=EyReal=EyImag=EzReal=EzImag=0.0;
     HxReal=HxImag=HyReal=HyImag=HzReal=HzImag=0.0;
     if (J)
      { 
        for(int nbf=0, no=0; no<NumObjects; no++)
         for(int nf=0; nf<Objects[no]->NumInteriorFaces; nf++, nbf++)
          { 
            cdouble EHScat[6];
            Get1BFFields(Objects[no], nf, Omega, X, EHScat);
            cdouble JAlpha = J->GetEntry(nbf);
            EHScat[0] *= JAlpha;
            EHScat[1] *= JAlpha;
            EHScat[2] *= JAlpha;
            EHScat[3] *= JAlpha;
            EHScat[4] *= JAlpha;
            EHScat[5] *= JAlpha;

            ExReal += real(EHScat[0]); ExImag += imag(EHScat[0]); 
            EyReal += real(EHScat[1]); EyImag += imag(EHScat[1]); 
            EzReal += real(EHScat[2]); EzImag += imag(EHScat[2]); 
            HxReal += real(EHScat[3]); HxImag += imag(EHScat[3]); 
            HyReal += real(EHScat[4]); HyImag += imag(EHScat[4]); 
            HzReal += real(EHScat[5]); HzImag += imag(EHScat[5]);

          }; // if (no...) if (nf...) 
      }; // if(J) ... 

     FMatrix->SetEntry(nr, 0, EHInc[0] + cdouble(ExReal, ExImag) );
     FMatrix->SetEntry(nr, 1, EHInc[1] + cdouble(EyReal, EyImag) );
     FMatrix->SetEntry(nr, 2, EHInc[2] + cdouble(EzReal, EzImag) );
     FMatrix->SetEntry(nr, 3, EHInc[3] + cdouble(HxReal, HxImag) );
     FMatrix->SetEntry(nr, 4, EHInc[4] + cdouble(HyReal, HyImag) );
     FMatrix->SetEntry(nr, 5, EHInc[5] + cdouble(HzReal, HzImag) );

   }; // for(int nr=0; nr<XMatrix->NR; nr++)

  return FMatrix;

}
  
/***************************************************************/
/* simpler interface to GetFields with only a single eval point*/
/***************************************************************/
void SWGGeometry::GetFields(IncField *IF, HVector *J, cdouble Omega, 
                            double *X, cdouble *EH)
{
  HMatrix XMatrix(1, 3, LHM_REAL, LHM_NORMAL, (void *)X);
  HMatrix FMatrix(1, 6, LHM_COMPLEX, LHM_NORMAL, (void *)EH);
  GetFields(IF, J, Omega, &XMatrix, &FMatrix);
} 

} // namespace buff
