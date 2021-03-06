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
 * TTaylorDuffy.cc:  implementation of the tetrahedron Taylor-Duffy method
 *                   for evaluating tetrahedra-product integrals over pairs
 *                   of tetrahedra with common vertices
 *
 * homer reid        1/2015
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <complex>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libTriInt.h>

#include "TTaylorDuffy.h"

using namespace scuff;

namespace scuff {
void GetScriptK(int WhichK, cdouble KParam, double X,
                int nMin, int nMax, cdouble KVector[7]);
                } 

#define DEFABSTOL 0.0
#define DEFRELTOL 1.0e-10
#define MAXFEVALS 10000
#define INTERVALS (MAXFEVALS/15)

#define II cdouble(0,1)

#define NUMREGIONS 18
#define NUMWPOWERS 8
#define MAXYDIM    5

#define MAXY1POW   4
#define MAXY2POW   4
#define MAXY3POW   4
#define MAXY4POW   2

#define NUMY1PS (MAXY1POW+1)
#define NUMY2PS (MAXY2POW+1)
#define NUMY3PS (MAXY3POW+1)
#define NUMY4PS (MAXY4POW+1)

typedef double UpsilonVector[NUMREGIONS][NUMWPOWERS][NUMY1PS][NUMY2PS][NUMY3PS][NUMY4PS];

typedef struct CoefficientRule
 { 
   unsigned char Region;
   unsigned char wPower;
   unsigned char y1Power;
   unsigned char y2Power;
   unsigned char y3Power;
   unsigned char y4Power;
   double Coefficient;

 } CoefficientRule;

/***************************************************************/
/* Data structure containing various data passed back and      */
/* forth among tetrahedron taylor-duffy routines.              */
/***************************************************************/
typedef struct TTDWorkspace
 {
   /* */
   int WhichCase;
   int NumPKs;
   int *PIndex;
   int *KIndex;
   cdouble *KParam;

   /* geometric data needed to compute X_d functions*/
   double RMatrix[6][6];

   /* coefficient rules */
   CoefficientRule *CRList[NUMPS];
   int NumCRs[NUMPS];
   int MinwPower[NUMPS], MaxwPower[NUMPS], MaxyPower[NUMPS][MAXYDIM];

   /* */
   bool NeedP[NUMPS];
   bool NeedK[NUMKS];

   int nCalls;

   FILE *IntegrandLogFile;

 } TTDWorkspace;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetXAndJPrime(TTDWorkspace *TTDW, const double *yVector,
                   double X[NUMREGIONS], double JPrime[NUMREGIONS]);

void GetScriptP(TTDWorkspace *TTDW, int WhichP,
                const double *yVector,
                double P[NUMREGIONS][NUMWPOWERS]);

void ComputeGeometricParameters(TTDArgStruct *Args,
                                TTDWorkspace *TTDW);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetUpsilon_CommonTet_TMuBDotBP_1(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTet_TMuBDotBP_2(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTet_TMuBDotBP_3(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTet_TMuBDotBP_4(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTet_TMuBDotBP_5(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTet_TMuBDotBP_6(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTriangle_TMuBDotBP_1(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTriangle_TMuBDotBP_2(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTriangle_TMuBDotBP_3(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTriangle_TMuBDotBP_4(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTriangle_TMuBDotBP_5(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonTriangle_TMuBDotBP_6(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonEdge_TMuBDotBP_1(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonEdge_TMuBDotBP_2(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonEdge_TMuBDotBP_3(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonEdge_TMuBDotBP_4(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonEdge_TMuBDotBP_5(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 
void GetUpsilon_CommonEdge_TMuBDotBP_6(TTDWorkspace *TTDW, int Mu, double L1[3], double L2[3], double L3[3], double L1P[3], double L2P[3], double L3P[3], double D[3], double DP[3], double V0mXTorque[3], UpsilonVector Upsilon); 

/***************************************************************/
/***************************************************************/
/***************************************************************/
int TTDIntegrand(unsigned ydim, const double *yVector, void *parms,
                 unsigned fdim, double *f)
{
  (void) fdim; // unused

  /*--------------------------------------------------------------*/
  /*- extract parameters from data structure ---------------------*/
  /*--------------------------------------------------------------*/
  TTDWorkspace *TTDW = (TTDWorkspace *)parms;
  int NumPKs        = TTDW->NumPKs;
  int *PIndex       = TTDW->PIndex;
  int *KIndex       = TTDW->KIndex;
  cdouble *KParam   = TTDW->KParam;
  TTDW->nCalls++;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nOffset = ydim;

  /*--------------------------------------------------------------*/
  /*- prefetch values of the X_d and J_d^\prime functions for all */
  /*- regions.                                                    */
  /*--------------------------------------------------------------*/
  double X[NUMREGIONS], JPrime[NUMREGIONS];
  GetXAndJPrime(TTDW, yVector, X, JPrime);

  /*--------------------------------------------------------------*/
  /*- prefetch values of the scriptP vector for all P functions  -*/
  /*- we will need                                               -*/
  /*--------------------------------------------------------------*/
  double P[NUMPS][NUMREGIONS][NUMWPOWERS];
  for(int np=0; np<NUMPS; np++)
   if (TTDW->NeedP[np])
    GetScriptP(TTDW, np, yVector, P[np]);

  /*--------------------------------------------------------------*/
  /*- assemble the integrand vector by adding all subregions and  */
  /*- all n-values                                                */
  /*--------------------------------------------------------------*/
  cdouble K[NUMREGIONS][NUMWPOWERS+5];
  cdouble *Sum=(cdouble *)f;
  for(int npk=0; npk<NumPKs; npk++)
   { 
     int np = PIndex[npk];
     int nMin = TTDW->MinwPower[ np ];
     int nMax = TTDW->MaxwPower[ np ];

     // attempt to skip the calculation of ScriptK
     // if we can use values from the previous loop iteration
     for(int d=0; d<NUMREGIONS; d++)
      GetScriptK( KIndex[npk], KParam[npk], X[d],
                  nMin+nOffset, nMax+nOffset, K[d]);

     Sum[npk]=0.0;
     for(int n=nMin; n<=nMax; n++)
      for(int d=0; d<NUMREGIONS; d++)
       Sum[npk] += JPrime[d] * P[np][d][n] * K[d][n+nOffset];

   };

  if (TTDW->IntegrandLogFile)
   { for(unsigned ny=0; ny<ydim; ny++)
      fprintf(TTDW->IntegrandLogFile,"%e ",yVector[ny]);
     fprintVecCR(TTDW->IntegrandLogFile,Sum,NumPKs);
   };
 
  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitTTDArgs(TTDArgStruct *Args)
{
  Args->Q=0;
  Args->QP=0;  
  Args->AbsTol=0.0;
  Args->RelTol=1.0e-6;
  Args->MaxEval=10000;
  Args->CubaturePoints=0;
  Args->IntegrandLogFile=0;
  Args->XTorque[0]=Args->XTorque[1]=Args->XTorque[2]=0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *CreateTTDWorkspace(TTDArgStruct *Args)
{
  TTDWorkspace *TTDW=(TTDWorkspace *)mallocEC(sizeof *TTDW);
  TTDW->WhichCase = Args->WhichCase;
  int NumPKs  = TTDW->NumPKs = Args->NumPKs;
  int *PIndex = TTDW->PIndex = Args->PIndex;
  int *KIndex = TTDW->KIndex = Args->KIndex;
                TTDW->KParam = Args->KParam;

  memset(TTDW->NeedP, 0, NUMPS*sizeof(bool));
  memset(TTDW->NeedK, 0, NUMKS*sizeof(bool));
  for(int npk=0; npk<NumPKs; npk++)
   { if ( PIndex[npk]<0 || PIndex[npk]>=NUMPS )
      ErrExit("invalid PIndex (%i) in TTaylorDuffy",PIndex[npk]);
     if ( KIndex[npk]<0 || KIndex[npk]>=NUMKS )
         ErrExit("invalid KIndex (%i) in TTaylorDuffy",KIndex[npk]);
     TTDW->NeedP[PIndex[npk]] = true;
     TTDW->NeedK[KIndex[npk]] = true;
   };

  memset(TTDW->NumCRs,0,NUMPS*sizeof(int));
  memset(TTDW->CRList,0,NUMPS*sizeof(TTDW->CRList[0]));

  TTDW->IntegrandLogFile = Args->IntegrandLogFile;

  return (void *)TTDW;

}

void DestroyTTDWorkspace(void *pTTDW)
{
  TTDWorkspace *TTDW=(TTDWorkspace *)pTTDW;
  if (TTDW==0) return;
 
  for(int np=0; np<NUMPS; np++)
   if (TTDW->CRList[np])
    free(TTDW->CRList[np]);
  free(TTDW);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TTaylorDuffy(TTDArgStruct *Args, void *pTTDW)
{
  /***************************************************************/
  /* unpack fields from argument structure ***********************/
  /***************************************************************/
  int WhichCase    = Args->WhichCase;

  int NumPKs         = Args->NumPKs;
  int *KIndex        = Args->KIndex;

  double AbsTol      = Args->AbsTol;
  double RelTol      = Args->RelTol;
  int MaxEval        = Args->MaxEval;
  int CubaturePoints = Args->CubaturePoints;

  /***************************************************************/
  /* initialize TTDW structure to pass data to integrand routines */
  /***************************************************************/
  bool OwnsTTDW=false;
  TTDWorkspace *TTDW=(TTDWorkspace *)pTTDW;
  if (TTDW==0)
   { OwnsTTDW=true;
     TTDW=(TTDWorkspace *)CreateTTDWorkspace(Args);
   };

  ComputeGeometricParameters(Args,TTDW);
  
  /***************************************************************/
  /* evaluate the 2-, 3-, 4-, or 5- dimensional cubature         */
  /***************************************************************/
  static double Lower[5]={0.0, 0.0, 0.0, 0.0, 0.0};
  static double Upper[5]={1.0, 1.0, 1.0, 1.0, 1.0};
  int fDim=2*NumPKs;
  double *dResult=(double *)(Args->Result);
  double *dError=(double *)(Args->Error);
  TTDW->nCalls=0;
  int IntegralDimension = 6 - WhichCase;

  if (IntegralDimension<=2 || CubaturePoints!=0)
   CCCubature(CubaturePoints, fDim, 
              TTDIntegrand, (void *)TTDW, IntegralDimension,
              Lower, Upper, MaxEval, AbsTol, RelTol,
              ERROR_INDIVIDUAL, dResult, dError);
  else
   hcubature(fDim, TTDIntegrand, (void *)TTDW, IntegralDimension,
             Lower, Upper, MaxEval, AbsTol, RelTol,
             ERROR_INDIVIDUAL, dResult, dError);

  Args->nCalls = TTDW->nCalls;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int npk=0; npk<NumPKs; npk++)
   if ( KIndex[npk]==TTD_HELMHOLTZ || KIndex[npk]==TTD_GRADHELMHOLTZ)
    Args->Result[npk] /= (4.0*M_PI);

  if (OwnsTTDW)
   DestroyTTDWorkspace( (void *) TTDW );
}

/***************************************************************/
/* uXi[0..6] = { u1, u2, u3, Xi3, Xi2, Xi1 }                   */
/***************************************************************/
void GetuXiAndJPrime(double w, const double *yVector,
                     int WhichCase,
                     double uXi[NUMREGIONS][6],
                     double JPrime[NUMREGIONS])
{ 
  double y1=yVector[0];
  double y2=yVector[1];
  double y3=0.0;
  double y4=0.0;

  if ( WhichCase<=TTD_COMMONTRIANGLE )
   y3=yVector[2]; 
  if ( WhichCase<=TTD_COMMONEDGE     )
   y4=yVector[3]; 

  double y1y2=y1*y2;
  double y1y3=y1*y3;
  double y1y2y3=y1*y2*y3;
  double y1y2y3y4=y1*y2*y3*y4;
  double y2y3=y2*y3;
  double y2y3y4=y2y3*y4;
  double y3y4=y3*y4;

  double wy1       = w*y1;
  double wy2       = w*y2;
  double wy1y2     = wy1*y2;
  double wy1y3     = wy1*y3;
  double wy1y2y3   = wy1y2*y3;
  double wy1y2y4   = wy1y2*y4;
  double wy1y2y3y4 = wy1y2y3*y4;

  double w1My1     = w*(1.0-y1);
  double wy11My2   = w*y1*(1.0-y2);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  switch(WhichCase) { 
     case TTD_COMMONTET:
      uXi[0][0]=wy1y2;     uXi[0][1]=wy1;        uXi[0][2]=w;
      uXi[1][0]=wy1;       uXi[1][1]=w;          uXi[1][2]=wy2;
      uXi[2][0]=wy1y2;     uXi[2][1]=wy1;        uXi[2][2]=-w1My1;
      uXi[3][0]=wy1;       uXi[3][1]=wy1y2;      uXi[3][2]=w*(1.0-y1+y1y2);
      uXi[4][0]=w;         uXi[4][1]=wy1;        uXi[4][2]=wy1y2;
      uXi[5][0]=wy1;       uXi[5][1]=wy1y2;      uXi[5][2]=-w1My1;
      uXi[6][0]=wy1y2;     uXi[6][1]=-wy11My2;   uXi[6][2]=w1My1;
      uXi[7][0]=w1My1;     uXi[7][1]=-wy1;       uXi[7][2]=-wy1y2;
      uXi[8][0]=wy11My2;   uXi[8][1]=-wy1y2;     uXi[8][2]=-w*(1.0-y1+y1y2);
      uXi[9][0]=-wy1y2;    uXi[9][1]=wy11My2;    uXi[9][2]=w*(1.0-y1y2);
      uXi[10][0]=-wy1y2;   uXi[10][1]=w-wy1y2;   uXi[10][2]=wy1*(1.0-y2);
      uXi[11][0]=-wy1y2;   uXi[11][1]=wy11My2;   uXi[11][2]=-w1My1;
      uXi[12][0]=-wy1;     uXi[12][1]=-wy1y2;    uXi[12][2]=w1My1;
      uXi[13][0]=-w;       uXi[13][1]=-wy1;      uXi[13][2]=-wy1y2;  
      uXi[14][0]=-wy1;     uXi[14][1]=-wy1y2;    uXi[14][2]=-w*(1.0-y1+y1y2);
      uXi[15][0]=-wy1y2;   uXi[15][1]=-wy1;      uXi[15][2]=w1My1;
      uXi[16][0]=-wy1;     uXi[16][1]=-w;        uXi[16][2]=-w*(1.0-y2);
      uXi[17][0]=-wy1y2;   uXi[17][1]=-wy1;      uXi[17][2]=-w;
      break;

     case TTD_COMMONTRIANGLE:
      uXi[0][0]=wy1y2y3;        uXi[0][1]=wy1y2;           uXi[0][2]=wy1;                 uXi[0][3]=w1My1;
      uXi[1][0]=wy1y2;          uXi[1][1]=wy1;             uXi[1][2]=wy1y3;               uXi[1][3]=w1My1;
      uXi[2][0]=wy1y2y3;        uXi[2][1]=wy1y2;           uXi[2][2]=-wy1*(1.0-y2);       uXi[2][3]=w*(1.0-y1y2);
      uXi[3][0]=wy1y2;          uXi[3][1]=wy1y2y3;         uXi[3][2]=wy1*(1.0-y2+y2y3);   uXi[3][3]=w1My1;
      uXi[4][0]=wy1;            uXi[4][1]=wy1y2;           uXi[4][2]=wy1y2y3;             uXi[4][3]=w1My1;
      uXi[5][0]=wy1y2;          uXi[5][1]=wy1y2y3;         uXi[5][2]=-wy11My2;            uXi[5][3]=w*(1.0-y1y2);
      uXi[6][0]=wy1y2y3;        uXi[6][1]=-wy1y2*(1.0-y3); uXi[6][2]=wy11My2;             uXi[6][3]=w1My1;
      uXi[7][0]=wy11My2;        uXi[7][1]=-wy1y2;          uXi[7][2]=-wy1y2y3;            uXi[7][3]=w*(1.0-y1+y1y2y3);
      uXi[8][0]=wy1y2*(1.0-y3); uXi[8][1]=-wy1y2y3;        uXi[8][2]=-wy1*(1.0-y2+y2y3);  uXi[8][3]=w*(1.0-y1y2+y1y2y3);
      uXi[9][0]=-wy1y2y3;       uXi[9][1]=wy1y2*(1.0-y3);  uXi[9][2]=wy1*(1.0-y2y3);      uXi[9][3]=w1My1;
      uXi[10][0]=-wy1y2y3;      uXi[10][1]=wy1*(1.0-y2y3); uXi[10][2]=wy1y2*(1.0-y3);     uXi[10][3]=w1My1;
      uXi[11][0]=-wy1y2y3;      uXi[11][1]=wy1y2*(1.0-y3); uXi[11][2]=-wy11My2;           uXi[11][3]=w*(1.0-y1y2);
      uXi[12][0]=-wy1y2;        uXi[12][1]=-wy1y2y3;       uXi[12][2]=wy11My2;            uXi[12][3]=w1My1;
      uXi[13][0]=-wy1;          uXi[13][1]=-wy1y2;         uXi[13][2]=-wy1y2y3;           uXi[13][3]=w*(1.0-y1+y1y2y3);
      uXi[14][0]=-wy1y2;        uXi[14][1]=-wy1y2y3;       uXi[14][2]=-wy1*(1.0-y2+y2y3); uXi[14][3]=w*(1.0-y1y2+y1y2y3);
      uXi[15][0]=-wy1y2y3;      uXi[15][1]=-wy1y2;         uXi[15][2]=wy1*(1.0-y2);       uXi[15][3]=w1My1;
      uXi[16][0]=-wy1y2;        uXi[16][1]=-wy1;           uXi[16][2]=-wy1*(1.0-y3);      uXi[16][3]=w*(1.0-y1y3);
      uXi[17][0]=-wy1y2y3;      uXi[17][1]=-wy1y2;         uXi[17][2]=-wy1;               uXi[17][3]=w;
      break;

     case TTD_COMMONEDGE:    
      uXi[0][0]=wy1y2y3y4;        uXi[0][1]=wy1y2y3;           uXi[0][2]=wy1y2;                 uXi[0][3]=wy11My2;                 uXi[0][4]=w-wy1y2y3;
      uXi[1][0]=wy1y2y4;          uXi[1][1]=wy1y2;             uXi[1][2]=wy1y2y3;               uXi[1][3]=wy11My2;                 uXi[1][4]=w-wy1y2;
      uXi[2][0]=wy1y2y3y4;        uXi[2][1]=wy1y2y3;           uXi[2][2]=-wy1y2*(1.0-y3);       uXi[2][3]=wy1*(1.0-y2y3);          uXi[2][4]=w-wy1y2y3;
      uXi[3][0]=wy1y2y3;          uXi[3][1]=wy1y2y3y4;         uXi[3][2]=wy1y2*(1.0-y3+y3y4);   uXi[3][3]=wy1*(1.0-y2);            uXi[3][4]=w-wy1y2y3;
      uXi[4][0]=wy1y2;            uXi[4][1]=wy1y2y3;           uXi[4][2]=wy1y2y3y4;             uXi[4][3]=wy11My2;                 uXi[4][4]=w-wy1y2;
      uXi[5][0]=wy1y2y3;          uXi[5][1]=wy1y2y3y4;         uXi[5][2]=-wy1y2*(1.0-y3);       uXi[5][3]=wy1*(1.0-y2y3);          uXi[5][4]=w-wy1y2y3;
      uXi[6][0]=wy1y2y3y4;        uXi[6][1]=-wy1y2y3*(1.0-y4); uXi[6][2]=wy1y2*(1.0-y3);        uXi[6][3]=wy11My2;                 uXi[6][4]=w-wy1y2y3y4;
      uXi[7][0]=wy1y2*(1.0-y3);   uXi[7][1]=-wy1y2y3;          uXi[7][2]=-wy1y2y3y4;            uXi[7][3]=wy1*(1.0-y2+y2y3y4);     uXi[7][4]=w-wy1y2+wy1y2y3;
      uXi[8][0]=wy1y2y3*(1.0-y4); uXi[8][1]=-wy1y2y3y4;        uXi[8][2]=-wy1y2*(1.0-y3+y3y4);  uXi[8][3]=wy1*(1.0-y2y3+y2y3y4);   uXi[8][4]=w*(1.0-y1y2y3+y1y2y3y4);
      uXi[9][0]=-wy1y2y3y4;       uXi[9][1]=wy1y2y3*(1.0-y4);  uXi[9][2]=wy1y2*(1.0-y3y4);      uXi[9][3]=wy11My2;                 uXi[9][4]=w-wy1y2y3;
      uXi[10][0]=-wy1y2y3y4;      uXi[10][1]=wy1y2*(1.0-y3y4); uXi[10][2]=wy1y2y3*(1.0-y4);     uXi[10][3]=wy11My2;                uXi[10][4]=w-wy1y2;
      uXi[11][0]=-wy1y2y3y4;      uXi[11][1]=wy1y2y3*(1.0-y4); uXi[11][2]=-wy1y2*(1.0-y3);      uXi[11][3]=wy1*(1.0-y2y3);         uXi[11][4]=w-wy1y2y3;
      uXi[12][0]=-wy1y2y3;        uXi[12][1]=-wy1y2y3y4;       uXi[12][2]=wy1y2*(1.0-y3);       uXi[12][3]=wy1*(1.0-y2);           uXi[12][4]=w*(1.0-y1y2y3+y1y2y3y4);
      uXi[13][0]=-wy1y2;          uXi[13][1]=-wy1y2y3;         uXi[13][2]=-wy1y2y3y4;           uXi[13][3]=wy1*(1.0-y2+y2y3y4);    uXi[13][4]=w*(1.0-y1y2+y1y2y3);
      uXi[14][0]=-wy1y2y3;        uXi[14][1]=-wy1y2y3y4;       uXi[14][2]=-wy1y2*(1.0-y3+y3y4); uXi[14][3]=wy1*(1.0-y2y3+y2y3y4);  uXi[14][4]=w*(1.0-y1y2y3+y1y2y3y4);
      uXi[15][0]=-wy1y2y3y4;      uXi[15][1]=-wy1y2y3;         uXi[15][2]=wy1y2*(1.0-y3);       uXi[15][3]=wy1*(1.0-y2);           uXi[15][4]=w;
      uXi[16][0]=-wy1y2y3;        uXi[16][1]=-wy1y2;           uXi[16][2]=-wy1y2*(1.0-y4);      uXi[16][3]=wy1*(1.0-y2*y4);        uXi[16][4]=w;
      uXi[17][0]=-wy1y2y3y4;      uXi[17][1]=-wy1y2y3;         uXi[17][2]=-wy1y2;               uXi[17][3]=wy1;                    uXi[17][4]=w;
      break;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  switch(WhichCase)
   { 
     case TTD_COMMONTET:
       JPrime[0]=y1;
       for(int d=1; d<NUMREGIONS; d++)
        JPrime[d]=JPrime[0];
       JPrime[1]=JPrime[16]=1.0;
       break;

     case TTD_COMMONTRIANGLE:
       JPrime[0]=y1*y1*y2;
       for(int d=1; d<NUMREGIONS; d++)
        JPrime[d]=JPrime[0];
       JPrime[1]=JPrime[16]=y1*y1;
       break;

     case TTD_COMMONEDGE:
       JPrime[0]=y1*y1*y1*y2*y2*y3;
       for(int d=1; d<NUMREGIONS; d++)
        JPrime[d]=JPrime[0];
       JPrime[1]=JPrime[16]=y1*y1*y1*y2*y2;
       break;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetXAndJPrime(TTDWorkspace *TTDW, const double *yVector,
                   double X[NUMREGIONS], double JPrime[NUMREGIONS])
{
  int WhichCase = TTDW->WhichCase;
  double uXi[NUMREGIONS][6];
  GetuXiAndJPrime(1.0, yVector, WhichCase, uXi, JPrime);
 
  int uXiDim = 7 - WhichCase;
  memset(X, 0, NUMREGIONS*sizeof(double));
  for(int m=0; m<uXiDim; m++)
   for(int n=m; n<uXiDim; n++)
    for(int d=0; d<NUMREGIONS; d++)
     X[d] += uXi[d][m] * TTDW->RMatrix[m][n] * uXi[d][n];

  for(int d=0; d<NUMREGIONS; d++)
   X[d] = sqrt(X[d]);
}

#if 0
#define THRESHOLD 1.0e-12
void GetnMinMax(double P[NUMREGIONS][NUMWPOWERS], int nMinMax[2])
{
  // 
  int nMin=100;
  int nMax=-1;
  for(int d=0; d<NUMREGIONS; d++)
   for(int n=0; n<NUMWPOWERS; n++)
    { if ( fabs(P[d][n]) > THRESHOLD )
       { if (n<nMin) nMin=n;
         if (n>nMax) nMax=n;
       };
    };
  nMinMax[0]=nMin;
  nMinMax[1]=nMax;

}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScriptP(TTDWorkspace *TTDW, int WhichP, const double *yVector,
                double P[NUMREGIONS][NUMWPOWERS])
{
  int WhichCase=TTDW->WhichCase;

  /***************************************************************/
  /* initialize vectors of powers of Y variables                 */
  /***************************************************************/
  double yPowers[MAXYDIM][NUMY1PS];
  yPowers[0][0]=yPowers[1][0]=yPowers[2][0]=yPowers[3][0]=1.0;
  int NumYs = 6 - WhichCase;
  for(int ny=0; ny<NumYs; ny++)
   for(int np=1; np<=TTDW->MaxyPower[WhichP][ny]; np++)
    yPowers[ny][np]=yPowers[ny][np-1]*yVector[ny];

  /***************************************************************/
  /* compute ScriptP polynomials for all subregions and all      */
  /* powers of w                                                 */
  /***************************************************************/
  for(int d=0; d<NUMREGIONS; d++)
   memset(P[d], 0, NUMWPOWERS*sizeof(double));

  for(int ncr=0; ncr<TTDW->NumCRs[WhichP]; ncr++)
   { 
     CoefficientRule *CR = &(TTDW->CRList[WhichP][ncr]);
     P[CR->Region][CR->wPower] += CR->Coefficient
                                  *yPowers[0][CR->y1Power]
                                  *yPowers[1][CR->y2Power]
                                  *yPowers[2][CR->y3Power]
                                  *yPowers[3][CR->y4Power];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ZeroUpsilonVector(UpsilonVector *Upsilon)
{
  for(int d=0; d<NUMREGIONS; d++)
   for(int w=0; w<NUMWPOWERS; w++)
    for(int y1p=0; y1p<NUMY1PS; y1p++)
     for(int y2p=0; y2p<NUMY2PS; y2p++)
      for(int y3p=0; y3p<NUMY3PS; y3p++)
       for(int y4p=0; y4p<NUMY4PS; y4p++)
        (*Upsilon)[d][w][y1p][y2p][y3p][y4p]=0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetCRList(UpsilonVector *Upsilon,
               int MinwPower, int MaxwPower, int MaxyPower[MAXYDIM],
               CoefficientRule **pCRList, int *pCRListLength)
{
#define UPSTHRESH 1.0e-10
  /***************************************************/
  /* first pass to count nonzero rules ***************/
  /***************************************************/
  int NumRules=0;
  for(int d=0; d<NUMREGIONS; d++)
   for(int wPower=MinwPower; wPower<=MaxwPower; wPower++)
    for(int y1p=0; y1p<=MaxyPower[0]; y1p++)
     for(int y2p=0; y2p<=MaxyPower[1]; y2p++)
      for(int y3p=0; y3p<=MaxyPower[2]; y3p++)
       for(int y4p=0; y4p<=MaxyPower[3]; y4p++)
        if ( UPSTHRESH < fabs( (*Upsilon)[d][wPower][y1p][y2p][y3p][y4p] ) )
         NumRules++;

  /***************************************************/
  /* second pass to populate list of rules ***********/
  /***************************************************/
  CoefficientRule *CRList = *pCRList;
  int CRListLength        = *pCRListLength;
  if ( CRList && CRListLength!=NumRules)
   { 
     Warn("%s:%i: mismatch in Taylor-Duffy workspace size (%i!=%i) (reallocating)...",__FILE__,__LINE__,CRListLength,NumRules);
     free(CRList);
     CRList=0;
   };

  if ( CRList==0 )
   { CRList=(CoefficientRule *)mallocEC(NumRules*sizeof(*CRList));
     *pCRList=CRList;
     *pCRListLength=NumRules;
   };
 
  for(int nr=0, d=0; d<NUMREGIONS; d++)
   for(int wPower=MinwPower; wPower<=MaxwPower; wPower++)
    for(int y1p=0; y1p<=MaxyPower[0]; y1p++)
     for(int y2p=0; y2p<=MaxyPower[1]; y2p++)
      for(int y3p=0; y3p<=MaxyPower[2]; y3p++)
       for(int y4p=0; y4p<=MaxyPower[3]; y4p++)
        if ( UPSTHRESH < fabs( (*Upsilon)[d][wPower][y1p][y2p][y3p][y4p] ) )
         { CRList[nr].Region=d;
           CRList[nr].wPower=wPower;
           CRList[nr].y1Power=y1p;
           CRList[nr].y2Power=y2p;
           CRList[nr].y3Power=y3p;
           CRList[nr].y4Power=y4p;
           CRList[nr].Coefficient=(*Upsilon)[d][wPower][y1p][y2p][y3p][y4p];
           nr++;
         };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeGeometricParameters(TTDArgStruct *Args,
                                TTDWorkspace *TTDW)
{ 
  double *V1       = Args->V1;
  double *V2       = Args->V2;
  double *V3       = Args->V3;
  double *V4       = Args->V4;
  double *V2P      = Args->V2P;
  double *V3P      = Args->V3P;
  double *V4P      = Args->V4P;
  double *Q        = Args->Q;
  double *QP       = Args->QP;
   
  /***************************************************************/
  /* the manual claims that certain parameters are not referenced*/
  /* in certain cases, which means the users might pass NULL for */
  /* those parameters, which could cause core dumps unless we do */
  /* the following                                               */
  /***************************************************************/
  if (Args->WhichCase==TTD_COMMONTET)
   { V2P=V2; V3P=V3; V4P=V4; } 
  else if (Args->WhichCase==TTD_COMMONTRIANGLE)
   { V2P=V2; V3P=V3; } 
  else if (Args->WhichCase==TTD_COMMONEDGE)
   { V2P=V2; }

  // R = x' - x 
  //   = u1*L1P + u2*L2P + u3*L3P + Xi1*DL1 + Xi2*DL2 + Xi3*DL3
  double L1[3], L2[3], L3[3], L1P[3], L2P[3], L3P[3];
  VecSub(V2,V1,L1);
  VecSub(V3,V2,L2);
  VecSub(V4,V3,L3);
  VecSub(V2P,V1,L1P);
  VecSub(V3P,V2P,L2P);
  VecSub(V4P,V3P,L3P);

  double DL1[3], DL2[3], DL3[3];
  VecSub(L1P,L1,DL1);
  VecSub(L2P,L2,DL2);
  VecSub(L3P,L3,DL3);

  TTDW->RMatrix[0][0] =     VecDot(L1P,L1P);
  TTDW->RMatrix[0][1] = 2.0*VecDot(L1P,L2P);
  TTDW->RMatrix[0][2] = 2.0*VecDot(L1P,L3P);
  TTDW->RMatrix[0][3] = 2.0*VecDot(L1P,DL3);
  TTDW->RMatrix[0][4] = 2.0*VecDot(L1P,DL2);
  TTDW->RMatrix[0][5] = 2.0*VecDot(L1P,DL1);

  TTDW->RMatrix[1][1] =     VecDot(L2P,L2P);
  TTDW->RMatrix[1][2] = 2.0*VecDot(L2P,L3P);
  TTDW->RMatrix[1][3] = 2.0*VecDot(L2P,DL3);
  TTDW->RMatrix[1][4] = 2.0*VecDot(L2P,DL2);
  TTDW->RMatrix[1][5] = 2.0*VecDot(L2P,DL1);

  TTDW->RMatrix[2][2] =     VecDot(L3P,L3P);
  TTDW->RMatrix[2][3] = 2.0*VecDot(L3P,DL3);
  TTDW->RMatrix[2][4] = 2.0*VecDot(L3P,DL2);
  TTDW->RMatrix[2][5] = 2.0*VecDot(L3P,DL1);

  TTDW->RMatrix[3][3] =     VecDot(DL3,DL3);
  TTDW->RMatrix[3][4] = 2.0*VecDot(DL3,DL2);
  TTDW->RMatrix[3][5] = 2.0*VecDot(DL3,DL1);

  TTDW->RMatrix[4][4] =     VecDot(DL2,DL2);
  TTDW->RMatrix[4][5] = 2.0*VecDot(DL2,DL1);

  TTDW->RMatrix[5][5] =     VecDot(DL1,DL1);

  double D[3]={0.0, 0.0, 0.0};
  double DP[3]={0.0, 0.0, 0.0};
  if ( Q && QP ) 
   { VecSub(V1,Q,D);
     VecSub(V1,QP,DP);
   };
  double L1dL1P = VecDot(L1, L1P);
  double L1dL2P = VecDot(L1, L2P);
  double L1dL3P = VecDot(L1, L3P);
  double L2dL1P = VecDot(L2, L1P);
  double L2dL2P = VecDot(L2, L2P);
  double L2dL3P = VecDot(L2, L3P);
  double L3dL1P = VecDot(L3, L1P);
  double L3dL2P = VecDot(L3, L2P);
  double L3dL3P = VecDot(L3, L3P);
  double L1dDP = VecDot(L1, DP);
  double L2dDP = VecDot(L2, DP);
  double L3dDP = VecDot(L3, DP);
  double DdL1P = VecDot(D, L1P);
  double DdL2P = VecDot(D, L2P);
  double DdL3P = VecDot(D, L3P);
  double DdDP  = VecDot(D, DP);

  double V0mXTorque[3];
  VecSub(V1, Args->XTorque, V0mXTorque);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(TTDW->CRList, 0, NUMPS*sizeof(CoefficientRule *));
  bool *NeedP   = TTDW->NeedP;
  int WhichCase = TTDW->WhichCase;
  UpsilonVector Upsilon;
  int MinwPower=0;
  int MaxwPower=0;
  int MaxyPower[MAXYDIM];
  if ( TTDW->NeedP[TTD_UNITY] )
   { 
     ZeroUpsilonVector(&Upsilon);

     if (WhichCase==TTD_COMMONTET) 
      { 
#include "UpsilonFiles/Upsilon_CommonTet_Unity.cc"
      }
     else if (WhichCase==TTD_COMMONTRIANGLE) 
      { 
#include "UpsilonFiles/Upsilon_CommonTriangle_Unity.cc" 
      }
     else if (WhichCase==TTD_COMMONEDGE)
      { 
#include "UpsilonFiles/Upsilon_CommonEdge_Unity.cc" 
      };

     int WhichP=TTD_UNITY;
     TTDW->MinwPower[WhichP]=MinwPower;
     TTDW->MaxwPower[WhichP]=MaxwPower;
     TTDW->MaxyPower[WhichP][0]=MaxyPower[0];
     TTDW->MaxyPower[WhichP][1]=MaxyPower[1];
     TTDW->MaxyPower[WhichP][2]=MaxyPower[2];
     TTDW->MaxyPower[WhichP][3]=MaxyPower[3];
     GetCRList(&Upsilon, MinwPower, MaxwPower, MaxyPower,
               TTDW->CRList + WhichP, TTDW->NumCRs + WhichP);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( TTDW->NeedP[TTD_BDOTBP] )
   { 
     ZeroUpsilonVector(&Upsilon);

     if (WhichCase==TTD_COMMONTET) 
      { 
#include "UpsilonFiles/Upsilon_CommonTet_BDotBP.cc"
      }
     else if (WhichCase==TTD_COMMONTRIANGLE) 
      { 
#include "UpsilonFiles/Upsilon_CommonTriangle_BDotBP.cc" 
      }
     else if (WhichCase==TTD_COMMONEDGE)
      { 
#include "UpsilonFiles/Upsilon_CommonEdge_BDotBP.cc" 
      };

     int WhichP=TTD_BDOTBP;
     TTDW->MinwPower[WhichP]=MinwPower;
     TTDW->MaxwPower[WhichP]=MaxwPower;
     TTDW->MaxyPower[WhichP][0]=MaxyPower[0];
     TTDW->MaxyPower[WhichP][1]=MaxyPower[1];
     TTDW->MaxyPower[WhichP][2]=MaxyPower[2];
     TTDW->MaxyPower[WhichP][3]=MaxyPower[3];
     GetCRList(&Upsilon, MinwPower, MaxwPower, MaxyPower,
               TTDW->CRList + WhichP, TTDW->NumCRs + WhichP);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#if 0
  for(int WhichP=TTD_RXUNITY; WhichP<=TTD_RZUNITY; WhichP++)
   if ( TTDW->NeedP[WhichP] )
    {  
       ZeroUpsilonVector(&Upsilon);

       int Mu=WhichP - TTD_RXUNITY; // Mu=0, 1, or 2 
       double L1Mu=L1[Mu], L2Mu=L2[Mu], L3Mu=L3[Mu];
       double L1PMu=L1P[Mu], L2PMu=L2P[Mu], L3PMu=L3P[Mu];

       if (WhichCase==TTD_COMMONTET) 
        { 
#include "UpsilonFiles/Upsilon_CommonTet_RMuUnity.cc"
        }
       else if (WhichCase==TTD_COMMONTRIANGLE) 
        { 
#include "UpsilonFiles/Upsilon_CommonTriangle_RMuUnity.cc"
        }
       else if (WhichCase==TTD_COMMONEDGE)
        { 
#include "UpsilonFiles/Upsilon_CommonEdge_RMuUnity.cc"
        };

       TTDW->MinwPower[WhichP]=MinwPower;
       TTDW->MaxwPower[WhichP]=MaxwPower;
       TTDW->MaxyPower[WhichP][0]=MaxyPower[0];
       TTDW->MaxyPower[WhichP][1]=MaxyPower[1];
       TTDW->MaxyPower[WhichP][2]=MaxyPower[2];
       TTDW->MaxyPower[WhichP][3]=MaxyPower[3];
       GetCRList(&Upsilon, MinwPower, MaxwPower, MaxyPower,
                 TTDW->CRList + WhichP, TTDW->NumCRs + WhichP);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int WhichP=TTD_RXBDOTBP; WhichP<=TTD_RZBDOTBP; WhichP++)
   if ( TTDW->NeedP[WhichP] )
    {  
       ZeroUpsilonVector(&Upsilon);

       int Mu=WhichP - TTD_RXBDOTBP; // Mu=0, 1, or 2 
       double L1Mu=L1[Mu], L2Mu=L2[Mu], L3Mu=L3[Mu];
       double L1PMu=L1P[Mu], L2PMu=L2P[Mu], L3PMu=L3P[Mu];

       if (WhichCase==TTD_COMMONTET) 
        { 
#include "UpsilonFiles/Upsilon_CommonTet_RMuBDotBP.cc"
        }
       else if (WhichCase==TTD_COMMONTRIANGLE) 
        { 
#include "UpsilonFiles/Upsilon_CommonTriangle_RMuBDotBP.cc"
        }
       else if (WhichCase==TTD_COMMONEDGE)
        { 
#include "UpsilonFiles/Upsilon_CommonEdge_RMuBDotBP.cc"
        };

       TTDW->MinwPower[WhichP]=MinwPower;
       TTDW->MaxwPower[WhichP]=MaxwPower;
       TTDW->MaxyPower[WhichP][0]=MaxyPower[0];
       TTDW->MaxyPower[WhichP][1]=MaxyPower[1];
       TTDW->MaxyPower[WhichP][2]=MaxyPower[2];
       TTDW->MaxyPower[WhichP][3]=MaxyPower[3];
       GetCRList(&Upsilon, MinwPower, MaxwPower, MaxyPower,
                 TTDW->CRList + WhichP, TTDW->NumCRs + WhichP);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int WhichP=TTD_TXUNITY; WhichP<=TTD_TZUNITY; WhichP++)
   if ( TTDW->NeedP[WhichP] )
    {  
       ZeroUpsilonVector(&Upsilon);

       int Mu=WhichP - TTD_TXUNITY; // Mu=0, 1, or 2 
       int MP1=(Mu+1)%3;
       int MP2=(Mu+2)%3;

       double L1MP1=L1[MP1], L2MP1=L2[MP1], L3MP1=L3[MP1];
       double L1MP2=L1[MP2], L2MP2=L2[MP2], L3MP2=L3[MP2];

       double L1PMP1=L1P[MP1], L2PMP1=L2P[MP1], L3PMP1=L3P[MP1];
       double L1PMP2=L1P[MP2], L2PMP2=L2P[MP2], L3PMP2=L3P[MP2];

       double V0mXTorqueMP1 = V0mXTorque[MP1];
       double V0mXTorqueMP2 = V0mXTorque[MP2];

     if (WhichCase==TTD_COMMONTET) 
      { 
#include "UpsilonFiles/Upsilon_CommonTet_TMuUnity.cc"
      }
     else if (WhichCase==TTD_COMMONTRIANGLE) 
      { 
#include "UpsilonFiles/Upsilon_CommonTriangle_TMuUnity.cc"
      }
     else if (WhichCase==TTD_COMMONEDGE)
      { 
#include "UpsilonFiles/Upsilon_CommonEdge_TMuUnity.cc"
      };

       TTDW->MinwPower[WhichP]=MinwPower;
       TTDW->MaxwPower[WhichP]=MaxwPower;
       TTDW->MaxyPower[WhichP][0]=MaxyPower[0];
       TTDW->MaxyPower[WhichP][1]=MaxyPower[1];
       TTDW->MaxyPower[WhichP][2]=MaxyPower[2];
       TTDW->MaxyPower[WhichP][3]=MaxyPower[3];
       GetCRList(&Upsilon, MinwPower, MaxwPower, MaxyPower,
                 TTDW->CRList + WhichP, TTDW->NumCRs + WhichP);

    }; //  if ( TTDW->NeedP[WhichP] )

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int WhichP=TTD_TXBDOTBP; WhichP<=TTD_TZBDOTBP; WhichP++)
   if ( TTDW->NeedP[WhichP] )
    {
       int Mu=WhichP-TTD_TXBDOTBP;
       ZeroUpsilonVector(&Upsilon);
       if (WhichCase==TTD_COMMONTET)
        { 
          MinwPower=0;
          MaxwPower=7;
          MaxyPower[0]=4;
          MaxyPower[1]=2;
          MaxyPower[2]=0;
          MaxyPower[3]=0;

          GetUpsilon_CommonTet_TMuBDotBP_1(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTet_TMuBDotBP_2(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTet_TMuBDotBP_3(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTet_TMuBDotBP_4(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTet_TMuBDotBP_5(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTet_TMuBDotBP_6(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
        }
       else if (WhichCase==TTD_COMMONTRIANGLE)
        { 
          MinwPower=0;
          MaxwPower=6;
          MaxyPower[0]=4;
          MaxyPower[1]=4;
          MaxyPower[2]=2;
          MaxyPower[3]=0;

          GetUpsilon_CommonTriangle_TMuBDotBP_1(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTriangle_TMuBDotBP_2(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTriangle_TMuBDotBP_3(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTriangle_TMuBDotBP_4(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTriangle_TMuBDotBP_5(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonTriangle_TMuBDotBP_6(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
        }
       else if (WhichCase==TTD_COMMONEDGE )
        { 
          MinwPower=0;
          MaxwPower=5;
          MaxyPower[0]=4;
          MaxyPower[1]=4;
          MaxyPower[2]=4;
          MaxyPower[3]=2;

          GetUpsilon_CommonEdge_TMuBDotBP_1(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonEdge_TMuBDotBP_2(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonEdge_TMuBDotBP_3(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonEdge_TMuBDotBP_4(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonEdge_TMuBDotBP_5(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
          GetUpsilon_CommonEdge_TMuBDotBP_6(TTDW, Mu, L1, L2, L3, L1P, L2P, L3P, D, DP, V0mXTorque, Upsilon);
        };

       TTDW->MinwPower[WhichP]=MinwPower;
       TTDW->MaxwPower[WhichP]=MaxwPower;
       TTDW->MaxyPower[WhichP][0]=MaxyPower[0];
       TTDW->MaxyPower[WhichP][1]=MaxyPower[1];
       TTDW->MaxyPower[WhichP][2]=MaxyPower[2];
       TTDW->MaxyPower[WhichP][3]=MaxyPower[3];
       GetCRList(&Upsilon, MinwPower, MaxwPower, MaxyPower,
                 TTDW->CRList + WhichP, TTDW->NumCRs + WhichP);

    };
#endif

} // GetGeometricParameters routine
