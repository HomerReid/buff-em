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
   bool *NeedDerivative;
   int NumDerivatives;
   double XTorque[3];
   int WhichKernel;
 } GMEData;
 
void GMEIntegrand(double *xA, double *bA, double DivbA,
                  double *xB, double *bB, double DivbB,
                  void *UserData, double *I)
{
  GMEData *Data        = (GMEData *)UserData;
  cdouble k            = Data->k;
  bool *NeedDerivative = Data->NeedDerivative;
  int NumDerivatives   = Data->NumDerivatives;
  int WhichKernel      = Data->WhichKernel;

  // fdim is the length of the output array in units of cdoubles
  int fdim;
  if (WhichKernel==KERNEL_STATIC)
   fdim = 3 + 4*NumDerivatives;
  else 
   fdim = 1 + NumDerivatives;

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
     Phi = exp(II*k*r) / (4.0*M_PI*r);
     PEFIE= DotProduct - ScalarProduct/(k*k);
     PEFIE = DotProduct - ScalarProduct/(k*k);
     zI[nf++] = PEFIE*Phi;
   }
  else if (WhichKernel==KERNEL_DESINGULARIZED)
   { 
     GetExpRelTable(IKR, 3, 4, ExpFac);
     Phi = ExpFac[0] / (4.0*M_PI*r);
     PEFIE = DotProduct - ScalarProduct/(k*k);
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

  if (!NeedDerivative)
   return;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double RFactor[6];
  RFactor[0]=R[0];
  RFactor[1]=R[1];
  RFactor[2]=R[2];

  // compute Tau = (X-XTorque) \cross R
  double XmXTorque[3];
  double *XTorque = Data->XTorque;
  XmXTorque[0] = xA[0] - XTorque[0];
  XmXTorque[1] = xA[1] - XTorque[1];
  XmXTorque[2] = xA[2] - XTorque[2];
  RFactor[3] = XmXTorque[1]*R[2]-XmXTorque[2]*R[1];
  RFactor[4] = XmXTorque[2]*R[0]-XmXTorque[0]*R[2];
  RFactor[5] = XmXTorque[0]*R[1]-XmXTorque[1]*R[0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (WhichKernel==KERNEL_STATIC)
   {
     double rPower[4];
     rPower[0] = -1.0/(4.0*M_PI*r*r*r);
     rPower[1] =  1.0/(8.0*M_PI*r);
     rPower[2] =  1.0/(12.0*M_PI);
     rPower[3] =    r/(24.0*M_PI);
     for(int Mu=0; Mu<6; Mu++)
      if (NeedDerivative[Mu])
       for(int nrPower=0; nrPower<4; nrPower++)
        {
          I[nf++] = RFactor[Mu]*DotProduct*rPower[nrPower];
          I[nf++] = RFactor[Mu]*ScalarProduct*rPower[nrPower];
        };
   }
  else 
   { cdouble Psi;
     if (WhichKernel==KERNEL_HELMHOLTZ)
      Psi=(IKR-1.0) * Phi / (r*r);
     else if (WhichKernel==KERNEL_DESINGULARIZED)
      Psi=(IKR-1.0)*ExpFac[1]/(4.0*M_PI*r*r*r);

     for(int Mu=0; Mu<6; Mu++)
      if (NeedDerivative[Mu])
       zI[nf++] = RFactor[Mu]*PEFIE*Psi;
   };

}

/***************************************************************/
/* Use the Taylor-Duffy method to compute the contribution of  */
/* a single tetrahedron-tetrahedron pair to the matrix         */
/* element of G and possibly its derivatives.                  */
/* GMETTI stands for 'G matrix element tet-tet integral.'      */
/***************************************************************/
void GetGMETTI_TaylorDuffy(SWGVolume *OA, int OVIA[4], int iQA,
                           SWGVolume *OB, int OVIB[4], int iQB,
                           int WhichKernel, bool *NeedDerivative,
                           cdouble k, int ncv, cdouble *TTI)
{
  /*-----------------------------------------------------------*/
  /*- derivatives vanish identically for the common-tet case  -*/
  /*-----------------------------------------------------------*/
  if (ncv==4 && NeedDerivative)
   { memset(TTI+1,0,6*sizeof(cdouble));
     NeedDerivative=0;
   };

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
  int PIndex[54], KIndex[54];
  cdouble KParam[54], TDI[54], Error[54];
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
     int npk=0;
     PIndex[npk] = TTD_BDOTBP; KIndex[npk] = TTD_HELMHOLTZ;  npk++;
     PIndex[npk] = TTD_UNITY;  KIndex[npk] = TTD_HELMHOLTZ;  npk++;
     for(int Mu=0; Mu<6; Mu++)
      if (NeedDerivative && NeedDerivative[Mu])
       { PIndex[npk] = TTD_RXBDOTBP + Mu;   KIndex[npk] = TTD_GRADHELMHOLTZ;  npk++;
         PIndex[npk] = TTD_RXUNITY  + Mu;   KIndex[npk] = TTD_GRADHELMHOLTZ;  npk++;
       };
     NumPKs=npk;
     for(int nnpk=0; nnpk<NumPKs; nnpk++)
      KParam[nnpk] = k;
   }
  else
   { 
     int npk=0;
     int rPowers[4]={-3,-1,0,1};
     for(int nrPower=1; nrPower<4; nrPower++)
      { PIndex[npk] = TTD_BDOTBP; KParam[npk] = rPowers[nrPower]; npk++;
        PIndex[npk] = TTD_UNITY;  KParam[npk] = rPowers[nrPower]; npk++;
      }
     for(int Mu=0; Mu<6; Mu++)
      if (NeedDerivative && NeedDerivative[Mu])
       for(int nrPower=0; nrPower<4; nrPower++)
        { PIndex[npk] = TTD_RXBDOTBP + Mu;  KParam[npk] = rPowers[nrPower]; npk++;
          PIndex[npk] = TTD_RXUNITY  + Mu;  KParam[npk] = rPowers[nrPower]; npk++;
        };
     NumPKs=npk;
     for(int nnpk=0; nnpk<NumPKs; nnpk++)
      KIndex[nnpk] = TTD_RP;
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
     TTI[0] = TDI[0] - NOK2*TDI[1];
     if (NeedDerivative)
      for(int Mu=0, nME=1, npk=2; Mu<6; Mu++)
       if (NeedDerivative[Mu])
        { TTI[nME++] = TDI[npk] - NOK2*TDI[npk+1];
          npk+=2;
        };
   }
  else
   { 
     double *dTTI = (double *)TTI;
     int npk=0;

     double GPreFactor[3]={1.0/(4.0*M_PI), 1.0/(4.0*M_PI), 1.0/(8.0*M_PI)};
     for(int nrPower=0; nrPower<3; nrPower++)
      { dTTI[npk] = GPreFactor[nrPower]*real(TDI[npk]);      npk++;
        dTTI[npk] = 9.0*GPreFactor[nrPower]*real(TDI[npk]);  npk++;
      };

     double dGPreFactor[4]={-1.0/(4.0*M_PI), 1.0/(8.0*M_PI), 1.0/(12.0*M_PI), 1.0/(24.0*M_PI)};
     for(int Mu=0; Mu<6; Mu++)
      if (NeedDerivative && NeedDerivative[Mu])
       for(int nrPower=0; nrPower<4; nrPower++)
        { dTTI[npk] = dGPreFactor[nrPower]*real(TDI[npk]);      npk++;
          dTTI[npk] = 9.0*dGPreFactor[nrPower]*real(TDI[npk]);  npk++;
        };
   };

}

/***************************************************************/
/* Get 'G' matrix elements via BF-BF integration, i.e.         */
/* by a single 6-dimensional integration over the full         */
/* supports of each basis function                             */
/***************************************************************/
void GetGME_BFBFInt(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB,
                    int WhichKernel, int NumDerivatives, 
                    bool *NeedDerivative,
                    cdouble Omega, cdouble *Result)
{
  GMEData MyData, *Data = &MyData;
  Data->k               = Omega;
  Data->WhichKernel     = WhichKernel;
  Data->NumDerivatives  = NumDerivatives;
  Data->NeedDerivative  = NeedDerivative;

  if (Data->NeedDerivative) 
   { Data->XTorque[0]=Data->XTorque[1]=Data->XTorque[2]=0.0;
     if (OA->OTGT) OA->OTGT->Apply(Data->XTorque);
     if (OA->GT)   OA->GT->Apply(Data->XTorque);
   };

  int fdim;
  if (WhichKernel==KERNEL_STATIC)
   fdim = 3 + 4*NumDerivatives;
  else 
   fdim = 1 + NumDerivatives;
  memset(Result, 0, 2*fdim*sizeof(double));

  int NumPts=4;
  cdouble Error[27];
  double Elapsed=Secs();
  BFBFInt(OA, nfA, OB, nfB, GMEIntegrand, (void *)Data, 2*fdim,
          (double *)Result, (double *)Error, NumPts, 0, 0.0);
  Elapsed=Secs()-Elapsed;
  //AddTaskTiming(0,Elapsed);
}

/***************************************************************/
/* Get 'G' matrix elements via Tetrahedron-Tetrahedron         */
/* integration, i.e. by 4 6-dimensional cubatures over pairs   */
/* of individual tetrahedra                                    */
/***************************************************************/
void GetGME_TetTetInt(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB,
                      int WhichKernel, int NumDerivatives, bool *NeedDerivative,
                      cdouble Omega, cdouble *Result)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GMEData MyData, *Data = &MyData;
  Data->k               = Omega;
  Data->WhichKernel     = WhichKernel;

  if (NumDerivatives>=0) 
   { Data->XTorque[0]=Data->XTorque[1]=Data->XTorque[2]=0.0;
     if (OA->OTGT) OA->OTGT->Apply(Data->XTorque);
     if (OA->GT)   OA->GT->Apply(Data->XTorque);
   };

  int fdim;
  if (WhichKernel==KERNEL_STATIC)
   fdim = 3 + 4*NumDerivatives;
  else 
   fdim = 1 + NumDerivatives;

  memset(Result, 0, 2*fdim*sizeof(double));

  /***************************************************************/
  /* do 4 6D cubatures to compute individual tet-tet contributions.*/
  /***************************************************************/
  SWGFace *FA = OA->Faces[nfA];
  SWGFace *FB = OB->Faces[nfB];
  cdouble TTI[27], Error[27];
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
         Data->NumDerivatives = ( (ncv==4) ? 0 : NumDerivatives );
         Data->NeedDerivative = ( (ncv==4) ? 0 : NeedDerivative );
         memset(TTI, 0, 27*sizeof(cdouble));
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
         memset(TTI,0,27*sizeof(cdouble));
         GetGMETTI_TaylorDuffy(OA, OVIA, iQA, OB, OVIB, iQB,
                               WhichKernel, NeedDerivative,
                               Omega, ncv, TTI);
         Elapsed=Secs()-Elapsed;
         //AddTaskTiming(ncv,Elapsed);
         for(int nf=0; nf<fdim; nf++)
          TTI[nf] *= 4.0*AreaFactor;
  
       }; // if (ncv<=1 ... else )
  
      double Sign = (ASign==BSign) ? 1.0 : -1.0;
      for(int nf=0; nf<fdim; nf++)
       Result[nf] += Sign*TTI[nf];

    }; // for(int ASign...for(int BSign...)

} // void SumTetTetInts(...)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeFIBBIData(SWGVolume *OA, int nfA,
                      SWGVolume *OB, int nfB,
                      double *GFI, double *dGFI)
{ 
  int NumDerivatives  = (dGFI ? 6 : 0);
  bool NeedDerivative[6];
  NeedDerivative[0] = (dGFI ? true : false);
  for(int Mu=1; Mu<6; Mu++)
   NeedDerivative[Mu] = NeedDerivative[0];

  cdouble zGFI[27];
  GetGME_TetTetInt(OA, nfA, OB, nfB, KERNEL_STATIC,
                   NumDerivatives, NeedDerivative, 0.0, zGFI);

  if (GFI)
   memcpy(GFI, (double *)zGFI, 6*sizeof(double));
  if (dGFI)
   { cdouble *zdGFI = zGFI + 3;
     memcpy(dGFI, (double *)zdGFI, 48*sizeof(double));
   };
}

/***************************************************************/
/* If dG is nonzero, it must point to an array of length 6.    */
/* The index Mu into this array has the following significance:*/
/*  Mu=0,1,2 -->  d/dx, d/dy, d/dz                             */
/*  Mu=3,4,5 -->  d/dTheta_x, d/dTheta_y, d/dTheta_z           */
/***************************************************************/
cdouble GetGMatrixElement(SWGVolume *OA, int nfA,
                          SWGVolume *OB, int nfB,
                          cdouble Omega, FIBBICache *GCache,
                          cdouble *dG, FIBBICache *dGCache)
{
  int NumDerivatives = (dG ? 6 : 0);
  bool NeedDerivative[6];
  NeedDerivative[0] = dG ? true : false;
  for(int Mu=1; Mu<6; Mu++)
   NeedDerivative[Mu]=NeedDerivative[0];

  /*--------------------------------------------------------------*/
  /* derivatives vanish identically for diagonal matrix elements  */
  /*--------------------------------------------------------------*/
  if ( dG && (OA==OB && nfA==nfB) )
   memset(dG, 0, 6*sizeof(cdouble));

  bool HaveCache = (GCache!=0);
  if (dG)
   HaveCache &= (dGCache!=0);

  /*--------------------------------------------------------------*/
  /*- count common vertices  -------------------------------------*/
  /*--------------------------------------------------------------*/
  double rRel;
  int ncv = CompareBFs(OA, nfA, OB, nfB, &rRel);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble GMEs[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  if ( ncv==0 )
   {
     GetGME_BFBFInt(OA, nfA, OB, nfB, KERNEL_HELMHOLTZ, 
                    NumDerivatives, NeedDerivative, Omega, GMEs);
   }
  else if (!HaveCache)
   { 
     GetGME_TetTetInt(OA, nfA, OB, nfB, KERNEL_HELMHOLTZ,
                      NumDerivatives, NeedDerivative, Omega, GMEs);

   }
  else
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     GetGME_BFBFInt(OA, nfA, OB, nfB, KERNEL_DESINGULARIZED,
                    NumDerivatives, NeedDerivative, Omega, GMEs);

     /***************************************************************/
     /* now look up or compute the desingularized contributions     */
     /* and add those in.                                           */
     /***************************************************************/
     double GFI[6];
     GCache->GetFIBBIData(OA, nfA, OB, nfB, GFI);
     cdouble k2=Omega*Omega;
     cdouble IK=II*Omega, IK2=IK*IK, IK3=IK2*IK, IK4=IK3*IK;
     GMEs[0] += (GFI[0]-GFI[1]/k2) + IK*(GFI[2]-GFI[3]/k2) + IK2*(GFI[4]-GFI[5]/k2);
     if (dG)
      { double dGFI[48];
        dGCache->GetFIBBIData(OA, nfA, OB, nfB, dGFI);
        for(int Mu=0; Mu<6; Mu++)
          { GMEs[Mu+1] +=        (dGFI[8*Mu+0] - dGFI[8*Mu+1]/k2)
                           + IK2*(dGFI[8*Mu+2] - dGFI[8*Mu+3]/k2)
                           + IK3*(dGFI[8*Mu+4] - dGFI[8*Mu+5]/k2)
                           + IK4*(dGFI[8*Mu+6] - dGFI[8*Mu+7]/k2);
          };
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (dG)
   memcpy(dG, GMEs+1, 6*sizeof(cdouble));
  return GMEs[0];

}

} // namespace buff
