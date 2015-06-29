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
  int nfA=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"nfA",                PA_INT,     1, 1, (void *)&nfA,            0, "nfA"},
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
  if (nfA==-1)
   nfA = lrand48() % O->NumInteriorFaces;
  printf("--nfA %i",nfA);

  double J[3];
  GetDQMoments(O, nfA, J, 0, false);
  double J2 = J[0]*J[0] + J[1]*J[1] + J[2]*J[2];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Omega=0.0;
  FILE *f=fopen("tPoynting.out","w");
  for( real(Omega)=1.0e-3; real(Omega)<=1.0; Omega*=exp(0.1*log(10.0)) )
   { 
     cdouble JGJ = GetGMatrixElement_SI(O, nfA, O, nfA, Omega, 0);
  
     double PRad = real(Omega) * J2 / (6.0*M_PI);

     fprintf(f,"%e %e %e %e \n",real(Omega),real(JGJ),imag(JGJ),PRad);
     fflush(f);
   };
  fclose(f);

}
