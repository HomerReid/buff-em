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
 * t1BFSITerms.cc
 *
 * homer reid       -- 6/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fenv.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;
using namespace buff;

#define II cdouble(0.0,1.0)

namespace buff {
void ExpRel23(cdouble x, cdouble *ExpRel2, cdouble *ExpRel3);
               }

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NFUN 3 
typedef struct IntegrandData
{ 
  double *XEval; 
  cdouble Omega;

} IntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void VIntegrand(double *XSource, double *b, double Divb,
                void *UserData, double *I)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  IntegrandData *Data = (IntegrandData *)UserData;
  double *XEval = Data->XEval;
  cdouble k     = Data->Omega;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double R[3], r;
  R[0] = XEval[0] - XSource[0];
  R[1] = XEval[1] - XSource[1];
  R[2] = XEval[2] - XSource[2];
  r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  cdouble ExpFac    = exp(II*k*r) / (4.0*M_PI*r);

/***************************************************************/
  double DeltaR;
  cdouble Gradient[3];
  cdouble ExpFacPD;
  for(int Mu=0; Mu<3; Mu++)
   { DeltaR = R[Mu]==0.0 ? 1.0e-4 : 1.0e-4*fabs(R[Mu]);
     R[Mu] += DeltaR;
     r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
     R[Mu] -= DeltaR;
     ExpFacPD = exp(II*k*r) / (4.0*M_PI*r);
     Gradient[Mu] = (ExpFacPD - ExpFac) / DeltaR;
   };
/***************************************************************/

  cdouble *zI = (cdouble *)I;
  zI[0] = Divb * Gradient[0] / (k*k);
  zI[1] = Divb * Gradient[1] / (k*k);
  zI[2] = Divb * Gradient[2] / (k*k);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SIntegrand(double *XSource, double *b, double Divb, double *nHat,
                void *UserData, double *I)
{
  IntegrandData *Data = (IntegrandData *)UserData;
  double *XEval = Data->XEval;
  cdouble k     = Data->Omega;

  double R[3];
  R[0] = XEval[0] - XSource[0];
  R[1] = XEval[1] - XSource[1];
  R[2] = XEval[2] - XSource[2];
  double RdnHat = R[0]*nHat[0] + R[1]*nHat[1] + R[2]*nHat[2]; 
  double r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble Xi = II*k*r;
  cdouble ER2, ER3;
  ExpRel23(Xi, &ER2, &ER3);
  cdouble h = ER2 / (4.0*M_PI*Xi);
  ExpRel23(-Xi, &ER2, &ER3);
  cdouble ExpFac = exp(Xi) / (4.0*M_PI*Xi);
  cdouble q = ExpFac * ER2 / (Xi*Xi);
  cdouble t = -1.0*ExpFac * (ER2 + 2.0*ER3) / (Xi*Xi*Xi*Xi);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble *zI = (cdouble *)I;
  zI[0] = Divb * (q*nHat[0] - k*k*RdnHat*t*R[0]) / (II*k);
  zI[1] = Divb * (q*nHat[1] - k*k*RdnHat*t*R[1]) / (II*k);
  zI[2] = Divb * (q*nHat[2] - k*k*RdnHat*t*R[2]) / (II*k);

}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{

  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble Omega=1.0;
  int nf=-1;
  int iQ=0;
  double DX[3]={2.0, 3.0, 4.0};
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nf",                 PA_INT,     1, 1, (void *)&nf,             0, "nfa"},
     {"iQ",                 PA_INT,     1, 1, (void *)&iQ,             0, "iQ"},
     {"DX",                 PA_DOUBLE,  3, 1, (void *)DX,              0, "DX"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry(GeoFile);
  SWGVolume *O   = G->Objects[0];
  
  srand48(time(0));
  if (nf==-1)
   nf = lrand48() % O->NumInteriorFaces;

  SWGFace *F     = O->Faces[nf];

  double XEval[3];
  XEval[0] = F->Centroid[0] + DX[0];
  XEval[1] = F->Centroid[1] + DX[1];
  XEval[2] = F->Centroid[2] + DX[2];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  IntegrandData MyData, *Data = &MyData;
  Data->XEval = XEval;
  Data->Omega = Omega;

  cdouble VResult[NFUN], VError[NFUN];
  TetInt(O, F->iPTet, iQ, 1.0, VIntegrand, Data, 2*NFUN,
         (double *)VResult, (double *)VError, 0, 0, 1.0e-4);

  cdouble SResult[NFUN], SError[NFUN];
  memset(SResult, 0, NFUN*sizeof(cdouble));
  memset(SError,  0, NFUN*sizeof(cdouble));
  for(int iF=0; iF<4; iF++)
   {
     cdouble PSResult[NFUN], PSError[NFUN];
     FaceInt(O, F->iPTet, iF, iQ, 1.0, SIntegrand, Data, 2*NFUN,
             (double *)PSResult, (double *)PSError,
             0, 0, 1.0e-4);
     for(int nf=0; nf<NFUN; nf++)
      SResult[nf] += PSResult[nf];
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nf=0; nf<NFUN; nf++)
   printf("%i:  %s  %s  %.1e\n",
           nf,CD2S(VResult[nf]),CD2S(SResult[nf]),
           RD(VResult[nf],SResult[nf]));

}
