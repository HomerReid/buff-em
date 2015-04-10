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
 * tTetTetInt
 *
 * homer reid       -- 6/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"
#include "TTaylorDuffy.h"
#include "MoreUtils.h"

using namespace scuff;
using namespace buff;

#define II cdouble(0.0,1.0)
#define NFUN 14

typedef struct IntegrandData
 { 
   int WhichKernel;
   cdouble k;
   double rPower;
   double xTorque[3];

 } IntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Integrand(double *xA, double *bA, double DivbA,
               double *xB, double *bB, double DivbB,
               void *UserData, double *I)
{
  IntegrandData *Data=(IntegrandData *)UserData;
  int WhichKernel = Data->WhichKernel;
  cdouble k       = Data->k;
  double rPower   = Data->rPower;
  double *xTorque = Data->xTorque;

  double R[3];
  R[0] = xA[0] - xB[0];
  R[1] = xA[1] - xB[1];
  R[2] = xA[2] - xB[2];
  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );

  double Tau[3];
  Tau[0] = (xA[1] - xTorque[1])*R[2] - (xA[2] - xTorque[2])*R[1];
  Tau[1] = (xA[2] - xTorque[2])*R[0] - (xA[0] - xTorque[0])*R[2];
  Tau[2] = (xA[0] - xTorque[0])*R[1] - (xA[1] - xTorque[1])*R[0];

  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB/9.0;

  cdouble K, Phi, Psi;
  if (WhichKernel==TTD_HELMHOLTZ)
   { cdouble ikr=II*k*r;
     Phi = exp(ikr) / (4.0*M_PI*r);
     Psi = (ikr-1.0) * Phi / (r*r);
   }
  else 
   { 
     Phi = Psi = pow(r,rPower);
   };

  cdouble *zI=(cdouble *)I;
  zI[0]  = Phi*DotProduct;
  zI[1]  = Phi*ScalarProduct;

  zI[2]  = Psi*R[0]*DotProduct;
  zI[3]  = Psi*R[0]*ScalarProduct;
  zI[4]  = Psi*R[1]*DotProduct;
  zI[5]  = Psi*R[1]*ScalarProduct;
  zI[6]  = Psi*R[2]*DotProduct;
  zI[7]  = Psi*R[2]*ScalarProduct;

  zI[8]  = Psi*Tau[0]*DotProduct;
  zI[9]  = Psi*Tau[0]*ScalarProduct;
  zI[10] = Psi*Tau[1]*DotProduct;
  zI[11] = Psi*Tau[1]*ScalarProduct;
  zI[12] = Psi*Tau[2]*DotProduct;
  zI[13] = Psi*Tau[2]*ScalarProduct;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetTetTetInt_TaylorDuffy(SWGVolume *O, int OVIA[4], int OVIB[4],
                              int ncv, double *QA, double *QB,
                              int WhichKernel, cdouble kParam,
                              int MaxEvals, double RelTol,
                              cdouble Result[14])
                              
{
  TTDArgStruct MyArgs, *Args=&MyArgs;
  InitTTDArgs(Args);

  Args->WhichCase=ncv;
  Args->V1     = O->Vertices             + 3*OVIA[0];
  Args->V2     = Args->V2P = O->Vertices + 3*OVIA[1];
  Args->V3     = Args->V3P = O->Vertices + 3*OVIA[2];
  Args->V4     = Args->V4P = O->Vertices + 3*OVIA[3];
  if (ncv<4)     Args->V4P = O->Vertices + 3*OVIB[3];
  if (ncv<3)     Args->V3P = O->Vertices + 3*OVIB[2];
  if (ncv<2)     Args->V2P = O->Vertices + 3*OVIB[1];
  Args->Q      = QA;
  Args->QP     = QB;
  if (O->OTGT) O->OTGT->Apply(Args->XTorque);
  if (O->GT) O->GT->Apply(Args->XTorque);

  int PIndex[NFUN], KIndex[NFUN];
  cdouble KParam[NFUN];

  Args->PIndex=PIndex;
  Args->KIndex=KIndex;
  Args->KParam=KParam;

  PIndex[0]  = TTD_BDOTBP;
  PIndex[1]  = TTD_UNITY;

  PIndex[2]  = TTD_RXBDOTBP;
  PIndex[3]  = TTD_RXUNITY;
  PIndex[4]  = TTD_RYBDOTBP;
  PIndex[5]  = TTD_RYUNITY;
  PIndex[6]  = TTD_RZBDOTBP;
  PIndex[7]  = TTD_RZUNITY;

  PIndex[8]  = TTD_TXBDOTBP;
  PIndex[9]  = TTD_TXUNITY;
  PIndex[10] = TTD_TYBDOTBP;
  PIndex[11] = TTD_TYUNITY;
  PIndex[12] = TTD_TZBDOTBP;
  PIndex[13] = TTD_TZUNITY;

  for(int npk=0; npk<14; npk++)
   { Args->KIndex[npk] = WhichKernel;
     Args->KParam[npk] = kParam;
   };
  Args->NumPKs=14;
  
  cdouble Error[14];
  Args->Result=Result;
  Args->Error=Error;
  Args->RelTol=RelTol;
  Args->MaxEval=MaxEvals;

  TTaylorDuffy(Args);

}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble Omega=1.0;
  int ntA=-1;
  int iQA=-1;
  int ntB=-1;
  int iQB=-1;
  int ncv=-1;
  cdouble k=0.0;
  double rPower=0.0;
  int BFOrder=0;
  int MaxBFEvals=100000;
  int MaxTDEvals=100000;
  double RelTolBF=1.0e-6;
  double RelTolTD=1.0e-6;
  
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFile,  0, ".buffgeo file"},
     {"ntA",        PA_INT,     1, 1, (void *)&ntA,      0, "ntA"},
     {"iQA",        PA_INT,     1, 1, (void *)&iQA,      0, "iQA"},
     {"ntB",        PA_INT,     1, 1, (void *)&ntB,      0, "ntB"},
     {"iQB",        PA_INT,     1, 1, (void *)&iQB,      0, "iQB"},
     {"ncv",        PA_INT,     1, 1, (void *)&ncv,      0, "ncv"},
     {"k",          PA_CDOUBLE, 1, 1, (void *)&k,        0, "k"},
     {"rPower",     PA_DOUBLE,  1, 1, (void *)&rPower,   0, "rPower"},
     {"BFOrder",    PA_INT,     1, 1, (void *)&BFOrder,    0, ""},
     {"MaxBFEvals", PA_INT,     1, 1, (void *)&MaxBFEvals, 0, ""},
     {"MaxTDEvals", PA_INT,     1, 1, (void *)&MaxTDEvals, 0, ""},
     {"RelTolBF",   PA_DOUBLE,  1, 1, (void *)&RelTolBF,   0, ""},
     {"RelTolTD",   PA_DOUBLE,  1, 1, (void *)&RelTolTD,   0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry(GeoFile);
  SWGVolume *O = G->Objects[0];
  
  srand48(time(0));
  int NT=O->NumTets;
  if (ntA==-1) ntA = lrand48() % O->NumTets;
  if (iQA==-1) iQA = lrand48() % 4;
  if (ntB==-1) ntB = lrand48() % O->NumTets;
  if (iQB==-1) iQB = lrand48() % 4;

  if (ncv!=-1) 
   { if (ncv==4) 
      ntB=ntA;
     else
      for (ntB=ntA+1; ntB!=ntA; ntB = (ntB+1)%NT)
       if (ncv==CompareTets(O,ntA,O,ntB))
        break;
   };

  int OVIA[4], OVIB[4];
  ncv=CompareTets(O,ntA,O,ntB,OVIA, OVIB);
  printf("--ntA %i --iQA %i ",ntA, iQA);
  printf("--ntB %i --iQB %i ",ntB, iQB);
  printf("--Omega %s \n",CD2S(Omega));
  printf("--ncv %i \n",ncv);
  
  int WhichKernel;
  cdouble KParam;
  if (k==0.0)
   { WhichKernel = TTD_RP;
     KParam = rPower;
   }
  else 
   { WhichKernel = TTD_HELMHOLTZ;
     KParam = k; 
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble ITD[NFUN];
  if( ncv > 1 )
   { double *QA = O->Vertices + 3*(O->Tets[ntA]->VI[iQA]);
     double *QB = O->Vertices + 3*(O->Tets[ntB]->VI[iQB]);
     GetTetTetInt_TaylorDuffy(O, OVIA, OVIB, ncv, QA, QB,
                              WhichKernel, KParam,
                              MaxTDEvals, RelTolTD, ITD);

     SWGFace *FA = O->Faces[ O->Tets[ntA]->FI[iQA] ];
     SWGFace *FB = O->Faces[ O->Tets[ntB]->FI[iQB] ];
     double AreaFactor = FA->Area * FB->Area;
     for(int nf=0; nf<NFUN; nf++)
      ITD[nf] *= 4.0*AreaFactor;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble IBF[NFUN], EBF[NFUN];
  IntegrandData MyData, *Data=&MyData;
  Data->WhichKernel = WhichKernel;
  Data->k = k;
  Data->rPower=rPower;
  memset(Data->xTorque,0,3*sizeof(double));
  if (O->OTGT) O->OTGT->Apply(Data->xTorque);
  if (O->GT)   O->GT  ->Apply(Data->xTorque);
  TetTetInt(O, ntA, iQA, 1.0, O, ntB, iQB, 1.0,
            Integrand, (void *)&Omega,
            2*NFUN, (double *)IBF, (double *)EBF,
            BFOrder, MaxBFEvals, RelTolBF);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Compare(IBF, ITD, 14, "BF", "TD");

}
