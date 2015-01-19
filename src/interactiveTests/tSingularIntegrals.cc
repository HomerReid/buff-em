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
 * tUMatrix.cc
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

using namespace scuff;
using namespace buff;

#define MAXSTR 1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace buff{
cdouble GetGMatrixElement_SI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, int Order);

              }

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GVIData 
 { cdouble k;
   double Eta;
 } GVIData;

void GVIntegrand(double *xA, double *bA, double DivbA,
                 double *xB, double *bB, double DivbB,
                 void *UserData, double *I)
{
  GVIData *Data = (GVIData *)UserData;

  cdouble k   = Data->k;
  double Eta = Data->Eta;

  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);

  double r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]) + Eta;

  cdouble ExpFac = exp(II*k*r)/(4.0*M_PI*r);

  double DotProduct = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];

  cdouble *zI = (cdouble *)I;
  zI[0] = (DotProduct - DivbA*DivbB/(k*k)) * ExpFac;
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
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nfa",                PA_INT,     1, 1, (void *)&nfa,            0, "nfa"},
     {"nfb",                PA_INT,     1, 1, (void *)&nfb,            0, "nfb"},
     {"ncv",                PA_INT,     1, 1, (void *)&ncv,            0, "ncv"},
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
  cdouble USI0 = GetGMatrixElement_SI(O, nfa, O, nfb, Omega, 0);
  cdouble USI20 = GetGMatrixElement_SI(O, nfa, O, nfb, Omega, 20);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetDefaultCD2SFormat("(%+.8e,%+.8e)");
  printf("USI0:     %s  \n",CD2S(USI0));
  printf("USI20:    %s \n",CD2S(USI20));
  GVIData MyData, *Data=&MyData;
  Data->k=Omega;
  double Eta;
  for(double Eta=0.01; Eta>=1.0e-10; Eta*=0.1)
   { Data->Eta = Eta;
     cdouble UVI, Error;
     BFBFInt(O, nfa, O, nfb, GVIntegrand, (void *)Data,
             2, (double *)&UVI, (double *)&Error, 0, 1000000, 1.0e-6);
     printf("UVI  :    %s (Eta=%e) (error=%s)\n",CD2S(UVI),Eta,CD2S(Error));
   };

}
