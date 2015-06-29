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
 * TetTetIntegral_SI.cc -- surface-integral (Bleszynski) computation
 *                      -- of tetrahedron-product integrals         
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

#define II cdouble(0.0,1.0)

/***************************************************************/
/* routine to compute the h, w, p functions described in the   */
/* memo                                                        */
/***************************************************************/
#define EXPRELTOL  1.0e-8
void Gethwp(cdouble x, cdouble hwp[3], cdouble *hwpPrime)
{
  if ( abs(x) >= 0.1 )
   { 
     cdouble x2=x*x;
     cdouble x3=x2*x;
     cdouble ExpRel1 = exp(x) - 1.0;
     cdouble ExpRel2 = ExpRel1 - x;
     cdouble ExpRel3 = ExpRel2 - 0.5*x2;

     hwp[0] = ExpRel2 / x;
     hwp[1] = ExpRel3;
     hwp[2] = hwp[0]/x  - hwp[1]/x3 - 1.0/3.0;

     if (hwpPrime)
      { cdouble x4=x3*x;
        hwpPrime[0] = ExpRel1/x - ExpRel2/x2;
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

     // in this loop, Term = x^n / n!
     cdouble Term = x2 / 2.0;
     for(double n=3.0; n<100.1; n+=1.0)
      { 
        Term *= x / n;

        hwp[0] += Term / (n+1.0);
        hwp[1] += Term;
        hwp[2] += Term / ( (n+3.0)*(n+1.0) );

        if (hwpPrime)
         { hwpPrime[0] += Term / (n+2.0);
           hwpPrime[2] += Term / ( (n+4.0)*(n+2.0) );
         };

        if ( abs(Term) < EXPRELTOL*abs(hwp[1]) )
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

void FFIntegrand_GMatrixElement(
                   double *xA, double *bA, double DivbA, double *nHatA,
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
  cdouble k2  = k*k;
  cdouble Xi  = II*k*r;

  cdouble hwp[3], hwpPrime[3];
  Gethwp(Xi, hwp, NeedDerivatives ? hwpPrime : 0 );

  cdouble h=hwp[0];
  cdouble w=hwp[1];
  cdouble p=hwp[2];

  double DotProd = (bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2]);
  cdouble PreFac  = -DivbA*DivbB / (9.0*k2);
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
        cdouble Factor = II*k*R[Mu]/r;

        T1 = Factor*DotProd*hP;
        T2 = Factor*PreFac*(9.0*hP + wP + k2*DQdotR*pP)
              +k2*PreFac*(QA[Mu]-QB[Mu])*p;
        zI[1+Mu] = NdotN * (T1 + T2) / (-4.0*M_PI*II*k);

        zI[4+Mu] = 0.0;
      };
   };

}

/***************************************************************/
/* Compute an integral over a tetrahedron product using the    */
/* surface-integral approach of Bleszynski et al.              */
/***************************************************************/
void TetTetInt_SI(SWGVolume *VA, int ntA, int iQa,
                  SWGVolume *VB, int ntB, int iQb,
                  UserFFIntegrand Integrand, void *UserData,
                  int fdim, double *Result, double *Error,
                  int Order, int MaxEvals, double RelTol)
{
  // 16 face-face integrals
  double *dResult = new double[fdim];
  double *dError = new double[fdim];
  memset(Result, 0, fdim*sizeof(double));
  memset(Error,  0, fdim*sizeof(double));
  for(int nfP=0; nfP<4; nfP++)
   for(int nfQ=0; nfQ<4; nfQ++)
    { 
      FaceFaceInt(VA, ntA, nfP, iQa, 1.0,
                  VB, ntB, nfQ, iQb, 1.0,
                  Integrand, UserData, fdim,
                  dResult, dError,
                  Order, MaxEvals, RelTol);

      for(int nf=0; nf<fdim; nf++)
       { Result[nf] += dResult[nf];
         Error[nf]  += dError[nf];
       };
    };
  delete[] dResult;
  delete[] dError;
}

} // namespace buff
