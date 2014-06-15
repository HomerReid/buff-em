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
 * t1BFFields.cc
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

namespace scuff{

void CalcGC(double R1[3], double R2[3],
            cdouble Omega, cdouble EpsR, cdouble MuR, 
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);

               }

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace buff{

void Get1BFFields_SI(SWGVolume *O, int nf, cdouble k, double X[3],
                     cdouble EH[6], int Order=0);

void Get1BFFields_DA(SWGVolume *O, int nf, cdouble k, double X[3],
                     cdouble EH[6], int NumTerms=1);

              }

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct VIData 
 { 
   double *x0;
   cdouble Omega;
 } VIData;

void VIntegrand(double *x, double *b, double Divb,
                 void *UserData, double *I)
{
  VIData *Data = (VIData *)UserData;

  cdouble G[3][3], C[3][3];
  CalcGC(Data->x0, x, Data->Omega, 1.0, 1.0, G, C, 0, 0);

  cdouble *zI = (cdouble *)I;
  cdouble EPrefac = II*(Data->Omega)*ZVAC;
  zI[0] = EPrefac * ( G[0][0]*b[0] + G[0][1]*b[1] + G[0][2]*b[2] );
  zI[1] = EPrefac * ( G[1][0]*b[0] + G[1][1]*b[1] + G[1][2]*b[2] );
  zI[2] = EPrefac * ( G[2][0]*b[0] + G[2][1]*b[1] + G[2][2]*b[2] );

  cdouble HPrefac = -II*(Data->Omega);
  zI[3] = HPrefac * ( C[0][0]*b[0] + C[0][1]*b[1] + C[0][2]*b[2] );
  zI[4] = HPrefac * ( C[1][0]*b[0] + C[1][1]*b[1] + C[1][2]*b[2] );
  zI[5] = HPrefac * ( C[2][0]*b[0] + C[2][1]*b[1] + C[2][2]*b[2] );
 
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
  int nf=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nf",                 PA_INT,     1, 1, (void *)&nf,             0, "nf"},
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
  if (nf==-1)
   nf = lrand48() % O->NumInteriorFaces;

  printf("--nf %i \n",nf);

  SWGFace *F = O->Faces[nf];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen("t1BFFields.dat","w");
  SetDefaultCD2SFormat("%.8e %.8e");
  for(double z=1.5*F->Radius; z<=100.0*F->Radius; z*=exp(0.05*log(10.0)) )
   { 
     double X[3];
     X[0] = F->Centroid[0];
     X[1] = F->Centroid[1];
     X[2] = F->Centroid[2] + z;

     cdouble EHVI[6], Error[6];
     VIData MyData, *Data=&MyData;
     Data->Omega=Omega;
     Data->x0=X;
     BFInt(O, nf, VIntegrand, (void *)Data,
           12, (double *)EHVI, (double *)Error, 0, 10000, 1.0e-6);
     double EVIAvg = abs(EHVI[0]) + abs(EHVI[1]) + abs(EHVI[2]);

     cdouble EHDI1[6], EHDI2[6];
     Get1BFFields_DA(O, nf, Omega, X, EHDI1, 1);
     double RDDI1 = (   abs( EHVI[0]-EHDI1[0] ) 
                      + abs( EHVI[1]-EHDI1[1] )
                      + abs( EHVI[2]-EHDI1[2] )  ) / EVIAvg;

     Get1BFFields_DA(O, nf, Omega, X, EHDI2, 2);
     double RDDI2 = (   abs( EHVI[0]-EHDI2[0] ) 
                      + abs( EHVI[1]-EHDI2[1] )
                      + abs( EHVI[2]-EHDI2[2] )  ) / EVIAvg;

     fprintf(f,"%e %e %e %e ",
                VecDistance(F->Centroid, X) / F->Radius,
                EVIAvg / 3.0, 
                RDDI1, RDDI2);
     fprintf(f,"%e %e %e %e %e %e ",
                real(EHVI[0]), imag(EHVI[0]),
                real(EHVI[1]), imag(EHVI[1]),
                real(EHVI[2]), imag(EHVI[2]));
     fprintf(f,"%e %e %e %e %e %e ",
                real(EHDI1[0]), imag(EHDI1[0]),
                real(EHDI1[1]), imag(EHDI1[1]),
                real(EHDI1[2]), imag(EHDI1[2]));
     fprintf(f,"%e %e %e %e %e %e ",
                real(EHDI2[0]), imag(EHDI2[0]),
                real(EHDI2[1]), imag(EHDI2[1]),
                real(EHDI2[2]), imag(EHDI2[2]));
     fprintf(f,"\n");


     fflush(f);
   };
  fclose(f);

}
