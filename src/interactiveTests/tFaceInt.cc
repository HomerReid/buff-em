/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * tFaceInt
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
void Integrand(double *x, double *b, double Divb, double *nHat,
               void *UserData, double *I)
{
  cdouble k = *(cdouble *)UserData;

  double r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );

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
  int nt=-1;
  int nf=-1;
  int iQ=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,  0, ".buffgeo file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,    0, "omega"},
     {"nt",                 PA_INT,     1, 1, (void *)&nt,       0, "nt"},
     {"nf",                 PA_INT,     1, 1, (void *)&nf,       0, "nf"},
     {"iQ",                 PA_INT,     1, 1, (void *)&iQ,       0, "iQ"},
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
  if (nt==-1)
   nt = lrand48() % O->NumTets;
  if (nf==-1)
   nf = lrand48() % 4;
  if (iQ==-1)
   iQ = lrand48() % O->NumVertices;
  printf("--nt %i --nf %i --iQ %i --omega %s\n",nt,nf,iQ,CD2S(Omega));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble IExact[NFUN], EExact[NFUN];
  FaceInt(O, nt, nf, iQ, 1.0, Integrand, &Omega,
          4, (double *)IExact, (double *)EExact,
          0, 0, 1.0e-8);

  int Orders[]={1, 4, 5, 7, 9, 13, 14, 16, 20, 25};
  #define NUMORDERS ( sizeof(Orders) / sizeof(Orders[0]) )
  cdouble IFixed[NUMORDERS][NFUN];

  for (int no=0; no<NUMORDERS; no++)
   FaceInt(O, nt, nf, iQ, 1.0, Integrand, &Omega,
           4, (double *)IFixed[no], 0, Orders[no], 0, 0.0);

  for (int nf=0; nf<NFUN; nf++)
   { printf("Function #%i: \n",nf);
     printf(" 0:  %s (%e) \n",CD2S(IExact[nf]),abs(EExact[nf]));
     for (int no=0; no<NUMORDERS; no++)
      printf("%2i: %s (%e) \n", Orders[no],
                                CD2S(IFixed[no][nf]),
                                RD(IFixed[no][nf], IExact[nf]));
   };

}
