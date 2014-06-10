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
                             cdouble Omega, int NumPts=0);

cdouble GetGMatrixElement_DA(SWGVolume *OA, int nfA, 
                             SWGVolume *OB, int nfB,
                             cdouble Omega, int NumTerms);
              }

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GVIntegrand(double *xA, double *bA, double DivbA,
                 double *xB, double *bB, double DivbB,
                 void *UserData, double *I)
{
  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);

  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);

  cdouble k = *((cdouble *)UserData);
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
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nfa",                PA_INT,     1, 1, (void *)&nfa,            0, "nfa"},
     {"nfb",                PA_INT,     1, 1, (void *)&nfb,            0, "nfb"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G1 = new SWGGeometry(GeoFile);
  SWGGeometry *G2 = new SWGGeometry(GeoFile);

  SWGVolume *O1 = G1->Objects[0];
  SWGVolume *O2 = G2->Objects[0];
  
  srand48(time(0));
  if (nfa==-1)
   nfa = lrand48() % O1->NumInteriorFaces;
  if (nfb==-1)
   nfb = lrand48() % O1->NumInteriorFaces;

  printf("--nfa %i --nfb %i\n",nfa,nfb);

  SWGFace *FA = O1->Faces[nfa];
  SWGFace *FB = O2->Faces[nfb];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen("tUElement.dat","w");
  SetDefaultCD2SFormat("%.8e %.8e");
  for(double z=0.1; z<=10.0; z*=exp(0.05*log(10.0)) )
   { 
     O2->Transform("DISPLACED 0 0 %e",z);

     cdouble UVI, Error;
     BFBFInt(O1, nfa, O2, nfb, GVIntegrand, (void *)&Omega, 
             2, (double *)&UVI, (double *)&Error, 0, 10000, 1.0e-6);

     cdouble USI   = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega);
     cdouble USI5  = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 5);
     cdouble USI14 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 14);
     cdouble USI20 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 20);
     cdouble USI25 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 25);
     cdouble UD1   = GetGMatrixElement_DA(O1, nfa, O2, nfb, Omega, 1);
     cdouble UD2   = GetGMatrixElement_DA(O1, nfa, O2, nfb, Omega, 2);

     printf("%.3e %.1e %.1e %.1e\n",z,abs(USI),abs(UVI),abs(USI)/abs(UVI));

     fprintf(f,"%e %e %e %e %e %e %e %e %e %e\n",
                VecDistance(FA->Centroid,FB->Centroid)
                 / fmax(FA->Radius, FB->Radius),
                real(USI), imag(USI),
                RD(USI, UVI),
                RD(USI, USI5),
                RD(USI, USI14),
                RD(USI, USI20),
                RD(USI, USI25),
                RD(USI, UD1), 
                RD(USI, UD2));
     fflush(f);

     O2->UnTransform();
   };
  fclose(f);

}
