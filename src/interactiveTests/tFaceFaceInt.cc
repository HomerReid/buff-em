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
 * tFaceFaceInt
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

#define II cdouble(0.0,1.0)
#define NFUN 2

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Integrand(double *xA, double *bA, double DivbA, double *nHatA,
               double *xB, double *bB, double DivbB, double *nHatB,
               void *UserData, double *I)
{
  cdouble k = *(cdouble *)UserData;

  double R[3];
  R[0] = xB[0] - xA[0];
  R[1] = xB[1] - xA[1];
  R[2] = xB[2] - xA[2];
  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );

  cdouble *zI=(cdouble *)I;
  zI[0] = 1.0;
  zI[1] = exp(II*k*r) / 4.0*M_PI*r;
  
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
  int nfA=-1;
  int iQA=-1;
  int ntB=-1;
  int nfB=-1;
  int iQB=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,  0, ".buffgeo file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,    0, "omega"},
     {"ntA",                PA_INT,     1, 1, (void *)&ntA,      0, "ntA"},
     {"nfA",                PA_INT,     1, 1, (void *)&nfA,      0, "nfA"},
     {"iQA",                PA_INT,     1, 1, (void *)&iQA,      0, "iQA"},
     {"ntB",                PA_INT,     1, 1, (void *)&ntB,      0, "ntB"},
     {"nfB",                PA_INT,     1, 1, (void *)&nfB,      0, "nfB"},
     {"iQB",                PA_INT,     1, 1, (void *)&iQB,      0, "iQB"},
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
  if (ntA==-1) ntA = lrand48() % O->NumTets;
  if (nfA==-1) nfA = lrand48() % 4; 
  if (iQA==-1) iQA = lrand48() % O->NumVertices;
  if (ntB==-1) ntB = lrand48() % O->NumTets;
  if (nfB==-1) nfB = lrand48() % 4; 
  if (iQB==-1) iQB = lrand48() % O->NumVertices;
  printf("--ntA %i --nfA %i --iQA %i ",ntA,nfA,iQA);
  printf("--ntB %i --nfB %i --iQB %i ",ntB,nfB,iQB);
  printf("--Omega %s\n",CD2S(Omega));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble IExact[NFUN], EExact[NFUN];
  FaceFaceInt(O, ntA, nfA, iQA, 1.0, O, ntB, nfB, iQB, 1.0,
              Integrand, &Omega,
              4, (double *)IExact, (double *)EExact,
              0, 0, 1.0e-8);

  int Orders[]={1, 4, 5, 7, 9, 13, 14, 16, 20, 25};
  #define NUMORDERS ( sizeof(Orders) / sizeof(Orders[0]) )
  cdouble IFixed[NUMORDERS][NFUN];

  for (int no=0; no<NUMORDERS; no++)
   FaceFaceInt(O, ntA, nfA, iQA, 1.0, O, ntB, nfB, iQB, 1.0,
               Integrand, &Omega,
               4, (double *)IFixed[no], 0,
               Orders[no], 0, 1.0e-8);

  for (int nf=0; nf<NFUN; nf++)
   { printf("Function #%i: \n",nf);
     printf(" 0:   %s (%e) \n",CD2S(IExact[nf]),abs(EExact[nf]));
     for (int no=0; no<NUMORDERS; no++)
      printf("%2i: %s (%e) \n", Orders[no],
                                CD2S(IFixed[no][nf]),
                                RD(IFixed[no][nf], IExact[nf]));
   };

}
