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

#define II cdouble(0.0,1.0)

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

void Get1BFFields_VI(SWGVolume *O, int nf, cdouble Omega, double X[3],
                     cdouble EH[6], int NumPts=0);

void Get1BFFields_SI(SWGVolume *O, int nf, cdouble k, double X[3],
                     cdouble EH[6], int Order=0);

void Get1BFFields_DA(SWGVolume *O, int nf, cdouble k, double X[3],
                     cdouble EH[6], int NumTerms=1);

              }
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
  double z=0.0;
  int Order = 20;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nf",                 PA_INT,     1, 1, (void *)&nf,             0, "nf"},
     {"z",                  PA_DOUBLE,  1, 1, (void *)&z,              0, "z"},
     {"Order",              PA_INT,     1, 1, (void *)&Order,          0, "Order"},
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
  if (z!=0.0)
   { 
     double X[3];
     X[0] = F->Centroid[0] + 3.0*z;
     X[1] = F->Centroid[1] + 2.0*z;
     X[2] = F->Centroid[2] + 1.0*z;

     cdouble EHVI[6];
     Get1BFFields_VI(O, nf, Omega, X, EHVI, 33);

     cdouble EHSI[6];
     Get1BFFields_SI(O, nf, Omega, X, EHSI, Order);
     printf("X: %s (%e) %s (%e) (%e) \n", CD2S(EHVI[0]), abs(EHVI[0]), CD2S(EHSI[0]), abs(EHSI[0]), abs(EHSI[0])/abs(EHVI[0]));
     printf("X: %s (%e) %s (%e) (%e) \n", CD2S(EHVI[1]), abs(EHVI[1]), CD2S(EHSI[1]), abs(EHSI[1]), abs(EHSI[1])/abs(EHVI[1]));
     printf("X: %s (%e) %s (%e) (%e) \n", CD2S(EHVI[2]), abs(EHVI[2]), CD2S(EHSI[2]), abs(EHSI[2]), abs(EHSI[2])/abs(EHVI[2]));
     exit(1);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen("t1BFFields.dat","w");
  SetDefaultCD2SFormat("%.8e %.8e");
  double Theta = M_PI*drand48();
  double Phi   = 2.0*M_PI*drand48();
  double SinTheta = sin(Theta);
  double CosTheta = cos(Theta);
  double SinPhi   = sin(Phi);
  double CosPhi   = cos(Phi);
  for(double rRel=0.01; rRel<=100.0; rRel*=exp(0.05*log(10.0)) )
   { 
     double X[3];
     X[0] = F->Centroid[0] + rRel*F->Radius*SinTheta*CosPhi;
     X[1] = F->Centroid[1] + rRel*F->Radius*SinTheta*SinPhi;
     X[2] = F->Centroid[2] + rRel*F->Radius*CosTheta;

     cdouble EHVI[6];
     Get1BFFields_VI(O, nf, Omega, X, EHVI, 33);

     cdouble EHSI[6];
     Get1BFFields_SI(O, nf, Omega, X, EHSI, 20);

     cdouble EHDA[6];
     Get1BFFields_DA(O, nf, Omega, X, EHDA, 2);

     fprintf(f,"%e ",rRel);
     fprintf(f,"%e %e ",real(EHVI[0]),imag(EHVI[0]));
     fprintf(f,"%e %e ",real(EHVI[1]),imag(EHVI[1]));
     fprintf(f,"%e %e ",real(EHVI[2]),imag(EHVI[2]));
     fprintf(f,"%e %e ",real(EHSI[0]),imag(EHSI[0]));
     fprintf(f,"%e %e ",real(EHSI[1]),imag(EHSI[1]));
     fprintf(f,"%e %e ",real(EHSI[2]),imag(EHSI[2]));
     fprintf(f,"%e %e ",real(EHDA[0]),imag(EHDA[0]));
     fprintf(f,"%e %e ",real(EHDA[1]),imag(EHDA[1]));
     fprintf(f,"%e %e ",real(EHDA[2]),imag(EHDA[2]));
     fprintf(f,"\n");
     fflush(f);
   };
  fclose(f);

}
