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
#include "libscuffInternals.h"
#include "FieldGrid.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define MAXFUNC 50

using namespace scuff;
namespace buff {

#define II cdouble(0,1)

/***************************************************************/
/* get the E and H fields due to a single SWG basis function   */
/* populated with unit strength                                */
/***************************************************************/
void Get1BFFields(SWGVolume *O, int nf, cdouble k, double X[3],
                  cdouble EH[6])
{
  SWGFace *F = O->Faces[nf];
  double rRel = VecDistance(X, F->Centroid) / F->Radius;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  
  GetDQMoments(B, nfB, JB, QB, true);
  
  
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
  else if ( (FMatrix->NR != XMatrix->NR) || (FMatrix->NC!=NumFuncs) ) 
   { Warn(" ** warning: wrong-size FMatrix passed to GetFields(); allocating new matrix");
     FMatrix=new HMatrix(XMatrix->NR, NumFuncs, LHM_COMPLEX);
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
                         reduction(+:ExReal, ExImag, EyReal, EyImag, ExReal, EzImag, \
                                     HxReal, HxImag, HyReal, HyImag, HxReal, HzImag)
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
         Get1BFFields(G->Objects[no], nf, k, X, EHScat);
         IF->GetFields(X, EHInc);

         cdouble JAlpha = J->GetEntry(nbf);
         EHTot[0] = JAlpha*EHScat[0] + EHInt[0];
         EHTot[1] = JAlpha*EHScat[1] + EHInt[1];
         EHTot[2] = JAlpha*EHScat[2] + EHInt[2];
         EHTot[3] = JAlpha*EHScat[3] + EHInt[3];
         EHTot[4] = JAlpha*EHScat[4] + EHInt[4];
         EHTot[5] = JAlpha*EHScat[5] + EHInt[5];

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
  GetFields(IF, J, Omega, &XMatrix, &FMatrix, 0);
} 

} // namespace scuff
