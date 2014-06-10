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
 * GetFields.cc  -- libscuff class methods for computing scattered
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

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define MAXFUNC 50

using namespace scuff;

namespace scuff{

void CalcGC(double R1[3], double R2[3],
            cdouble Omega, cdouble EpsR, cdouble MuR, 
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

               }

cdouble ExpRel(cdouble x, int n);
namespace buff {

#define II cdouble(0,1)

/***************************************************************/
/* get 1BF fields using surface-integral method ****************/
/***************************************************************/
void Get1BFFields_SI(SWGVolume *O, int nf, cdouble k, double X[3],
                     cdouble EH[6], int Order=0)
{
}

/***************************************************************/
/* get 1BF fields using dipole / quadrupole approximation      */
/***************************************************************/
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

/***************************************************************/
/* get the E and H fields due to a single SWG basis function   */
/* populated with unit strength                                */
/***************************************************************/
void Get1BFFields(SWGVolume *O, int nf, cdouble k, double X[3],
                  cdouble EH[6])
{
  SWGFace *F = O->Faces[nf];
  double rRel = VecDistance(X, F->Centroid) / F->Radius;

  if (rRel > 20.0)
   Get1BFFields_DA(O, nf, k, X, EH, 1);
  else if (rRel > 10.0)
   Get1BFFields_DA(O, nf, k, X, EH, 2);
  else
   Get1BFFields_SI(O, nf, k, X, EH, 20);
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
     ExReal=ExImag=EyReal=EyImag=EzReal=EzImag=0.0;
     HxReal=HxImag=HyReal=HyImag=HzReal=HzImag=0.0;
     for(int no=0, nbf=0; no<NumObjects; no++)
      for(int nf=0; nf<Objects[no]->NumInteriorFaces; nf++, nbf++)
       { 
         cdouble EHScat[6], EHInc[6], EHTot[6];
         Get1BFFields(Objects[no], nf, Omega, X, EHScat);
         IF->GetFields(X, EHInc);

         cdouble JAlpha = J->GetEntry(nbf);
         EHTot[0] = JAlpha*EHScat[0] + EHInc[0];
         EHTot[1] = JAlpha*EHScat[1] + EHInc[1];
         EHTot[2] = JAlpha*EHScat[2] + EHInc[2];
         EHTot[3] = JAlpha*EHScat[3] + EHInc[3];
         EHTot[4] = JAlpha*EHScat[4] + EHInc[4];
         EHTot[5] = JAlpha*EHScat[5] + EHInc[5];

         ExReal += real(EHTot[0]); ExImag += imag(EHTot[0]); 
         EyReal += real(EHTot[1]); EyImag += imag(EHTot[1]); 
         EzReal += real(EHTot[2]); EzImag += imag(EHTot[2]); 
         HxReal += real(EHTot[3]); HxImag += imag(EHTot[3]); 
         HyReal += real(EHTot[4]); HyImag += imag(EHTot[4]); 
         HzReal += real(EHTot[5]); HzImag += imag(EHTot[5]);

       }; // for (int no=) ... for (int nf=...)

     FMatrix->SetEntry(nr, 0, cdouble(ExReal, ExImag) );
     FMatrix->SetEntry(nr, 1, cdouble(EyReal, EyImag) );
     FMatrix->SetEntry(nr, 2, cdouble(EzReal, EzImag) );
     FMatrix->SetEntry(nr, 3, cdouble(HxReal, HxImag) );
     FMatrix->SetEntry(nr, 4, cdouble(HyReal, HyImag) );
     FMatrix->SetEntry(nr, 5, cdouble(HzReal, HzImag) );
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

} // namespace scuff
