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
 * tGMatrixElement.cc
 *
 * homer reid       -- 2/2015
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"
#include "MoreUtils.h"

using namespace scuff;
using namespace buff;

#define MAXSTR 1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GBFData 
 { cdouble k;
   int rPower;
 } GBFData;

void GBFIntegrand(double *xA, double *bA, double DivbA,
                  double *xB, double *bB, double DivbB,
                  void *UserData, double *I)
{
  GBFData *Data = (GBFData *)UserData;

  cdouble k  = Data->k;
  int rPower = Data->rPower;

  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);
  
  double r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);

  cdouble Phi, Psi;
  if (rPower!=-1)
   { 
     Phi=Psi=pow(r,rPower);
   }
  else
   { 
     Phi = exp(II*k*r)/(4.0*M_PI*r);
     Psi = (II*k - 1.0) * Phi / (r*r);
   };

  double DotProduct = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  cdouble PEFIE = DotProduct - DivbA*DivbB/(k*k);

  double *XmXT=xA; // assume torque center is the origin

  cdouble *zI = (cdouble *)I;
  zI[0] = PEFIE * Phi;
  zI[1] = R[0] * PEFIE * Psi;
  zI[2] = R[1] * PEFIE * Psi;
  zI[3] = R[2] * PEFIE * Psi;
  zI[4] = (XmXT[1]*R[2] - XmXT[2]*R[1]) * PEFIE * Psi;
  zI[5] = (XmXT[2]*R[0] - XmXT[0]*R[2]) * PEFIE * Psi;
  zI[6] = (XmXT[0]*R[1] - XmXT[1]*R[0]) * PEFIE * Psi;

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
  int nfa=-1;
  int nfb=-1;
  int ncv=-1;
  int rPower=0;
  int MaxBFEvals=100000;
  double RelTol=1.0e-6;
  bool ForceBF=false;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nfa",                PA_INT,     1, 1, (void *)&nfa,            0, "nfa"},
     {"nfb",                PA_INT,     1, 1, (void *)&nfb,            0, "nfb"},
     {"ncv",                PA_INT,     1, 1, (void *)&ncv,            0, "ncv"},
     {"rPower",             PA_INT,     1, 1, (void *)&rPower,         0, "rPower"},
     {"MaxBFEvals",         PA_INT,     1, 1, (void *)&MaxBFEvals,     0, "max integrand evaluations for BF integration"},
     {"RelTol",             PA_DOUBLE,  1, 1, (void *)&RelTol,         0, "relative tolerance for BF integration"},
     {"ForceBF",            PA_BOOL,    0, 1, (void *)&ForceBF,        0, "force BF"},
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
  if (nfa==-1)
   nfa = lrand48() % O->NumInteriorFaces;
  if (nfb==-1)
   nfb = lrand48() % O->NumInteriorFaces;

  if (ncv!=1)
   { for(nfb=0; nfb<O->NumInteriorFaces-1; nfb++)
      if (CompareBFs(O,nfa,O,nfb,0)==ncv) 
       break; 
   };

  ncv=CompareBFs(O,nfa,O,nfb,0);
  printf("--nfa %i --nfb %i ",nfa,nfb); 
  printf(" (%i common vertices)\n",ncv);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble GHR[7];
  bool NeedDerivatives[6]={true, true, true, true, true, true};
  GHR[0]=GetGMatrixElement(O, nfa, O, nfb, Omega, NeedDerivatives,
                           GHR+1, rPower, ForceBF);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GBFData MyData, *Data=&MyData;
  Data->k=Omega;
  Data->rPower=rPower;
  cdouble GBF[7], Error[7];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  TetTetInt(O, O->Faces[nfa]->iPTet, O->Faces[nfa]->PIndex, 1.0, 
            O, O->Faces[nfb]->iPTet, O->Faces[nfb]->PIndex, 1.0, 
            GBFIntegrand, (void *)Data,
            14, (double *)GBF, (double *)Error, 0, MaxBFEvals, RelTol);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*
  BFBFInt(O, nfa, O, nfb, GBFIntegrand, (void *)Data,
          14, (double *)GBF, (double *)Error, 0, MaxBFEvals, RelTol);
*/

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Compare(GBF, GHR, 7, "BF", "HR");

}
