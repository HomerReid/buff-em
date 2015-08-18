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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

/*
 * GMatrixElements.cc -- libSWG routines for computing matrix elements
 *                    -- of the dyadic Helmholtz operator and its derivatives
 *                    -- between SWG basis functions
 *
 * homer reid         -- 5/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"
#include "TTaylorDuffy.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

using namespace scuff;
namespace buff {

#define II cdouble(0.0,1.0)

#define KERNEL_HELMHOLTZ 0
#define KERNEL_DESINGULARIZED 1
#define KERNEL_STATIC    2

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL  
void GetExpRelTable(cdouble x, int nMin, int nMax, cdouble *ExpRelTable)
{
  // Note: the comments below assume nMin=4, nMax=7 just for concreteness

  int n;
  cdouble Term, Sum;

  // compute (but do not sum) the 0th, 1st, ... 3rd terms in the 
  // series expansion of exp(x)
  for(Term=1.0, n=1; n<nMin; n++)
   Term*=x/((double)n);
  
  // set ExpRelTable[0] = x^4/4!
  //     ExpRelTable[1] = x^5/25
  // ...
  //     ExpRelTable[3] = x^7/7!
  for(n=nMin; n<=nMax; n++)
   { Term*=x/((double)n);
     ExpRelTable[n-nMin]=Term;
   };

  // compute Sum = \sum_{n=8}^\infty x^n/n!
  for(Sum=0.0 ; n<100; n++)
   { Term*=x/((double)n);
     Sum+=Term;
     double mag2Term=norm(Term), mag2Sum=norm(Sum);
     if ( mag2Term < EXPRELTOL2*mag2Sum )
      break;
   };
  
  // 
  ExpRelTable[nMax-nMin]+=Sum;
  for(n=nMax-1; n>=nMin; n--)
   ExpRelTable[n-nMin] += ExpRelTable[n-nMin+1];
} 

/***************************************************************/
/* integrand routine and user data structure for computing     */
/* G matrix elements using BFBFInt or TetTetInt                */
/***************************************************************/
typedef struct GMEData
 {
   cdouble k;
   int WhichKernel;
 } GMEData;
 
void GMEIntegrand(double *xA, double *bA, double DivbA,
                  double *xB, double *bB, double DivbB,
                  void *UserData, double *I)
{
  GMEData *Data    = (GMEData *)UserData;
  cdouble k        = Data->k;
  int WhichKernel  = Data->WhichKernel;
  int fdim         = (WhichKernel==KERNEL_STATIC) ? 3 : 1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);
  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );
  if ( r < 1.0e-12 )
   { memset(I, 0, 2*fdim*sizeof(double));
     return;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int nf=0;
  cdouble PEFIE=0.0, Phi=0.0, ExpFac[2], *zI = (cdouble *)I;

  cdouble IKR = II*k*r;
  if (WhichKernel==KERNEL_HELMHOLTZ)
   { 
     Phi      = exp(II*k*r) / (4.0*M_PI*r);
     PEFIE    = DotProduct - ScalarProduct/(k*k);
     zI[nf++] = PEFIE*Phi;
   }
  else if (WhichKernel==KERNEL_DESINGULARIZED)
   { 
     GetExpRelTable(IKR, 3, 3, ExpFac);
     Phi      = ExpFac[0] / (4.0*M_PI*r);
     PEFIE    = DotProduct - ScalarProduct/(k*k);
     zI[nf++] = PEFIE*Phi;
   }
  else // (WhichKernel==KERNEL_STATIC)
   {
     double rPower[3];
     rPower[0] = 1.0/(4.0*M_PI*r);
     rPower[1] = 1.0/(4.0*M_PI);
     rPower[2] =   r/(8.0*M_PI);
     I[nf++] = rPower[0]*DotProduct;
     I[nf++] = rPower[0]*ScalarProduct;
     I[nf++] = rPower[1]*DotProduct;
     I[nf++] = rPower[1]*ScalarProduct;
     I[nf++] = rPower[2]*DotProduct;
     I[nf++] = rPower[2]*ScalarProduct;
   };

}

/***************************************************************/
/* Use the Taylor-Duffy method to compute the contribution of  */
/* a single tetrahedron-tetrahedron pair to the matrix         */
/* element of G.                                               */
/***************************************************************/
cdouble GetGMETTI_TaylorDuffy(SWGVolume *OA, int OVIA[4], int iQA,
                              SWGVolume *OB, int OVIB[4], int iQB,
                              int WhichKernel, 
                              cdouble k, int ncv, cdouble *TTI=0)
{
  /*-----------------------------------------------------------*/
  /* initialize Taylor-Duffy argument structure                */
  /*-----------------------------------------------------------*/
  TTDArgStruct MyArgs, *Args=&MyArgs;
  InitTTDArgs(Args);

  /***************************************************************/
  /* specify tetrahedra geometry *********************************/
  /***************************************************************/
  Args->WhichCase=ncv;
  Args->V1     = OA->Vertices             + 3*OVIA[0];
  Args->V2     = Args->V2P = OA->Vertices + 3*OVIA[1];
  Args->V3     = Args->V3P = OA->Vertices + 3*OVIA[2];
  Args->V4     = Args->V4P = OA->Vertices + 3*OVIA[3];
  if (ncv<4)     Args->V4P = OB->Vertices + 3*OVIB[3];
  if (ncv<3)     Args->V3P = OB->Vertices + 3*OVIB[2];
  if (ncv<2)     Args->V2P = OB->Vertices + 3*OVIB[1];
  Args->Q      = OA->Vertices + 3*iQA;
  Args->QP     = OB->Vertices + 3*iQB;
  if (OA->OTGT) OA->OTGT->Apply(Args->XTorque);
  if (OA->GT) OA->GT->Apply(Args->XTorque);
 
  /***************************************************************/
  /* specify the components of the Taylor-Duffy integrand vector */
  /***************************************************************/
  int PIndex[6], KIndex[6];
  cdouble KParam[6], TDI[6], Error[6];
  int NumPKs;
  Args->PIndex  = PIndex;
  Args->KIndex  = KIndex;
  Args->KParam  = KParam;
  Args->Result  = TDI;
  Args->Error   = Error;
  Args->RelTol  = SWGGeometry::TaylorDuffyTolerance;
  Args->MaxEval = SWGGeometry::MaxTaylorDuffyEvals;

  if (WhichKernel==KERNEL_HELMHOLTZ)
   { 
     PIndex[0] = TTD_BDOTBP; KIndex[0] = TTD_HELMHOLTZ; KParam[0]=k;
     PIndex[1] = TTD_UNITY;  KIndex[1] = TTD_HELMHOLTZ; KParam[1]=k;
     NumPKs=2;
   }
  else // static kernels for FIBBIs
   { 
     PIndex[0] = TTD_BDOTBP; KIndex[0] = TTD_RP; KParam[0]=-1.0;
     PIndex[1] = TTD_UNITY;  KIndex[1] = TTD_RP; KParam[1]=-1.0;
     PIndex[2] = TTD_BDOTBP; KIndex[2] = TTD_RP; KParam[2]= 0.0;
     PIndex[3] = TTD_UNITY;  KIndex[3] = TTD_RP; KParam[3]= 0.0;
     PIndex[4] = TTD_BDOTBP; KIndex[4] = TTD_RP; KParam[4]= 1.0;
     PIndex[5] = TTD_UNITY;  KIndex[5] = TTD_RP; KParam[5]= 1.0;
     NumPKs=6;
   };

  /***************************************************************/
  /* calculate taylor-duffy integrals                            */
  /***************************************************************/
  Args->NumPKs  = NumPKs;
  TTaylorDuffy(Args);

  /***************************************************************/
  /* assemble results into output quantities                     */
  /***************************************************************/
  if (WhichKernel==KERNEL_HELMHOLTZ)
   { 
     cdouble NOK2=9.0/(k*k);
     return TDI[0] - NOK2*TDI[1];
   }
  else
   {
     if (TTI==0) return 0.0;

     double *dTTI = (double *)TTI;
     int npk=0;

     double GPreFactor[3]={1.0/(4.0*M_PI), 1.0/(4.0*M_PI), 1.0/(8.0*M_PI)};
     for(int nrPower=0; nrPower<3; nrPower++)
      { dTTI[npk] = GPreFactor[nrPower]*real(TDI[npk]);      npk++;
        dTTI[npk] = 9.0*GPreFactor[nrPower]*real(TDI[npk]);  npk++;
      };

     return 0.0;

   };

}

/***************************************************************/
/* Get 'G' matrix elements via BF-BF integration, i.e.         */
/* by a single 6-dimensional integration over the full         */
/* supports of each basis function                             */
/***************************************************************/
cdouble GetGME_BFBFInt(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB,
                       int WhichKernel, cdouble Omega, cdouble *Result=0)
{
  GMEData MyData, *Data = &MyData;
  Data->k               = Omega;
  Data->WhichKernel     = WhichKernel;

  int fdim = (WhichKernel==KERNEL_STATIC) ? 3 : 1;
  cdouble ResultBuffer[3];
  if (Result==0) Result=ResultBuffer;
  memset(Result, 0, 2*fdim*sizeof(double));

  int NumPts=4;
  cdouble Error[3];
  BFBFInt(OA, nfA, OB, nfB, GMEIntegrand, (void *)Data, 2*fdim,
          (double *)Result, (double *)Error, NumPts, 0, 0.0);

  return Result[0];
}

/***************************************************************/
/* Get 'G' matrix elements via Tetrahedron-Tetrahedron         */
/* integration, i.e. by 4 6-dimensional cubatures over pairs   */
/* of individual tetrahedra                                    */
/***************************************************************/
cdouble GetGME_TetTetInt(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB,
                      int WhichKernel, cdouble Omega, cdouble *Result=0)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GMEData MyData, *Data = &MyData;
  Data->k               = Omega;
  Data->WhichKernel     = WhichKernel;

  int fdim = (WhichKernel==KERNEL_STATIC) ? 3 : 1;
  cdouble ResultBuffer[3];
  if (Result==0) Result=ResultBuffer;
  memset(Result, 0, 2*fdim*sizeof(double));

  /***************************************************************/
  /* do 4 6D cubatures to compute individual tet-tet contributions.*/
  /***************************************************************/
  SWGFace *FA = OA->Faces[nfA];
  SWGFace *FB = OB->Faces[nfB];
  cdouble TTI[3], Error[3];
  double AreaFactor=FA->Area * FB->Area;
  for(int ASign=0; ASign<2; ASign++)
   for(int BSign=0; BSign<2; BSign++)
    { 
      int ntA = (ASign==0) ? FA->iPTet : FA->iMTet;
      int ntB = (BSign==0) ? FB->iPTet : FB->iMTet;
      int OVIA[4], OVIB[4];
      int ncv=CompareTets(OA, ntA, OB, ntB, OVIA, OVIB);
  
      /***************************************************************/
      /* compute the tet--tet integral using brute-force cubature or */
      /* taylor-duffy                                                */
      /***************************************************************/
      if (ncv<=1 || WhichKernel==KERNEL_DESINGULARIZED)
       { 
         int NumPts=16;
         int iQA = (ASign==0) ? FA->PIndex : FA->MIndex;
         int iQB = (BSign==0) ? FB->PIndex : FB->MIndex;
         double Elapsed=Secs();
         memset(TTI, 0, fdim*sizeof(cdouble));
         TetTetInt(OA, ntA, iQA, 1.0, OB, ntB, iQB, 1.0,
                   GMEIntegrand, (void *)Data, 2*fdim,
                   (double *)TTI, (double *)Error, NumPts, 0, 0);
         Elapsed=Secs()-Elapsed;
         //AddTaskTiming(1,Elapsed);
       }
      else
       { 
         int iQA = (ASign==0) ? FA->iQP : FA->iQM;
         int iQB = (BSign==0) ? FB->iQP : FB->iQM;
         double Elapsed=Secs();
         memset(TTI,0,fdim*sizeof(cdouble));
         GetGMETTI_TaylorDuffy(OA, OVIA, iQA, OB, OVIB, iQB,
                               WhichKernel, Omega, ncv, TTI);
         Elapsed=Secs()-Elapsed;
         //AddTaskTiming(ncv,Elapsed);
         for(int nf=0; nf<fdim; nf++)
          TTI[nf] *= 4.0*AreaFactor;
  
       }; // if (ncv<=1 ... else )
  
      double Sign = (ASign==BSign) ? 1.0 : -1.0;
      for(int nf=0; nf<fdim; nf++)
       Result[nf] += Sign*TTI[nf];

    }; // for(int ASign...for(int BSign...)

  return Result[0];

} // void SumTetTetInts(...)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeFIBBIData(SWGVolume *OA, int nfA,
                      SWGVolume *OB, int nfB,
                      double GFI[FIBBIDATALEN])
{ 
  GetGME_TetTetInt(OA, nfA, OB, nfB, KERNEL_STATIC, 0.0, (cdouble *)GFI);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetGMatrixElement(SWGVolume *OA, int nfA,
                          SWGVolume *OB, int nfB,
                          cdouble Omega, FIBBICache *GCache)
{
  bool HaveCache = (GCache!=0);

  /*--------------------------------------------------------------*/
  /*- count common vertices  -------------------------------------*/
  /*--------------------------------------------------------------*/
  double rRel;
  int ncv = CompareBFs(OA, nfA, OB, nfB, &rRel);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble GME=0.0;
  if ( ncv==0 )
   {
     GME=GetGME_BFBFInt(OA, nfA, OB, nfB, KERNEL_HELMHOLTZ, Omega);
   }
  else if (!HaveCache)
   { 
     GME=GetGME_TetTetInt(OA, nfA, OB, nfB, KERNEL_HELMHOLTZ, Omega);
   }
  else
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     GME=GetGME_BFBFInt(OA, nfA, OB, nfB, KERNEL_DESINGULARIZED, Omega);

     /***************************************************************/
     /* now look up or compute the desingularized contributions     */
     /* and add those in.                                           */
     /***************************************************************/
     double GFI[6];
     GCache->GetFIBBIData(OA, nfA, OB, nfB, GFI);
     cdouble k2=Omega*Omega;
     cdouble IK=II*Omega, IK2=IK*IK;
     GME += (GFI[0]-GFI[1]/k2) + IK*(GFI[2]-GFI[3]/k2) + IK2*(GFI[4]-GFI[5]/k2);
   };

  return GME;

}

} // namespace buff
