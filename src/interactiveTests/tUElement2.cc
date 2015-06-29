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
 * tUElement.cc
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


/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace buff{
cdouble GetGMatrixElement_SI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, int NumPoints);

cdouble GetGMatrixElement_DA(SWGVolume *OA, int nfA, 
                             SWGVolume *OB, int nfB,
                             cdouble Omega, int NumTerms);
              }

#define II cdouble(0.0,1.0)

int CountCommonVertices(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB)
{ 
  int ncv=0;
  SWGFace *FA = OA->Faces[nfA];
  SWGFace *FB = OB->Faces[nfB];

  if (FA->iQP == FB->iQP ) ncv++;
  if (FA->iQP == FB->iV1 ) ncv++;
  if (FA->iQP == FB->iV2 ) ncv++;
  if (FA->iQP == FB->iV3 ) ncv++;
  if (FA->iQP == FB->iQM ) ncv++;

  if (FA->iV1 == FB->iQP ) ncv++;
  if (FA->iV1 == FB->iV1 ) ncv++;
  if (FA->iV1 == FB->iV2 ) ncv++;
  if (FA->iV1 == FB->iV3 ) ncv++;
  if (FA->iV1 == FB->iQM ) ncv++;

  if (FA->iV2 == FB->iQP ) ncv++;
  if (FA->iV2 == FB->iV1 ) ncv++;
  if (FA->iV2 == FB->iV2 ) ncv++;
  if (FA->iV2 == FB->iV3 ) ncv++;
  if (FA->iV2 == FB->iQM ) ncv++;

  if (FA->iV3 == FB->iQP ) ncv++;
  if (FA->iV3 == FB->iV1 ) ncv++;
  if (FA->iV3 == FB->iV2 ) ncv++;
  if (FA->iV3 == FB->iV3 ) ncv++;
  if (FA->iV3 == FB->iQM ) ncv++;

  if (FA->iQM == FB->iQP ) ncv++;
  if (FA->iQM == FB->iV1 ) ncv++;
  if (FA->iQM == FB->iV2 ) ncv++;
  if (FA->iQM == FB->iV3 ) ncv++;
  if (FA->iQM == FB->iQM ) ncv++;

  return ncv;

}

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

  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble Omega=1.0;
  int nfa=-1;
  int nfb=-1;
  int ncv=-1;
  double rRel=0.0;
  double DX[3]={0.0, 0.0, 0.0};
  bool DoVI=false;
  bool RefSI=false;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nfa",                PA_INT,     1, 1, (void *)&nfa,            0, "nfa"},
     {"nfb",                PA_INT,     1, 1, (void *)&nfb,            0, "nfb"},
     {"ncv",                PA_INT,     1, 1, (void *)&ncv,            0, "ncv"},
     {"rRel",               PA_DOUBLE,  1, 1, (void *)&rRel,           0, "rRel"},
     {"DX",                 PA_DOUBLE,  3, 1, (void *)DX,              0, "DX"},
     {"DoVI",               PA_BOOL,    0, 1, (void *)&DoVI,           0, "do the volume integral"},
     {"RefSI",              PA_BOOL,    0, 1, (void *)&RefSI,          0, "use the SI calculation as the reference"},
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

  if (ncv!=-1)
   { if (ncv<0 || ncv>5) 
      ErrExit("invalid ncv value");
     for(nfb=0; nfb<O1->NumInteriorFaces; nfb++)
      if( CountCommonVertices(O1, nfa, O1, nfb)==ncv )
       break;
     if (nfb==O1->NumInteriorFaces)
      ErrExit("couldn't find a BF pair with %i common vertices",ncv);
   };

  if (rRel!=0.0)
   { 
     for(nfb=0; nfb<O1->NumInteriorFaces; nfb++)
      { double MyrRel = VecDistance(O1->Faces[nfa]->Centroid, O1->Faces[nfb]->Centroid) 
                       / fmax(O1->Faces[nfa]->Radius, O1->Faces[nfb]->Radius);
        if ( 0.75*rRel< MyrRel && MyrRel < 1.25*rRel )
         break;
      };
     if (nfb==O1->NumInteriorFaces)
      ErrExit("couldn't find a BF pair with rRel=%e",rRel);
   };

  O2->Transform("DISPLACED %e %e %e",DX[0],DX[1],DX[2]);
  
  ncv=CountCommonVertices(O1, nfa, O1, nfb);
  rRel=VecDistance(O1->Faces[nfa]->Centroid, O2->Faces[nfb]->Centroid)
        / fmax(O1->Faces[nfa]->Radius, O1->Faces[nfb]->Radius);
  printf("\n");
  printf("--nfa %i --nfb %i \n",nfa,nfb);
  printf("%i common vertices\n",ncv);
  printf("rRel = %e \n",rRel);
  printf("\n");

  FILE *f=fopen("tUElement2.pp","w");
  O1->DrawBF(nfa, f);
  O2->DrawBF(nfb, f);
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/

  cdouble UVI=0.0, Error;
  if (DoVI)
   BFBFInt(O1, nfa, O2, nfb, GVIntegrand, (void *)&Omega, 
           2, (double *)&UVI, (double *)&Error, 0, 10000, 1.0e-8);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble UVI16=0.0;
  double TimeVI16 = 0.0;
  if (DoVI)
   {
     Tic();
     for(int n=0; n<5; n++)
      BFBFInt(O1, nfa, O2, nfb, GVIntegrand, (void *)&Omega, 
              2, (double *)&UVI16, (double *)&Error, 16, 10000, 1.0e-8);
     TimeVI16 = Toc() / 5; 
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Tic();
  cdouble USI;
  for(int n=0; n<1; n++)
   USI = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 0);
  double Time0 = Toc() / 1; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble USI5;
  Tic();
  for(int n=0; n<100; n++)
   USI5 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 5);
  double Time5 = Toc() / 100; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble USI9;
  Tic();
  for(int n=0; n<20; n++)
   USI9 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 9);
  double Time9 = Toc() / 20; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble USI16;
  Tic();
  for(int n=0; n<20; n++)
   USI16 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 16);
  double Time16 = Toc() / 20; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble USI20;
  Tic();
  for(int n=0; n<1; n++)
   USI20 = GetGMatrixElement_SI(O1, nfa, O2, nfb, Omega, 20);
  double Time20 = Toc() / 1; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Tic();
  cdouble UD1;
  for(int n=0; n<100; n++)
   UD1 = GetGMatrixElement_DA(O1, nfa, O2, nfb, Omega, 1); 
  double TimeD1 = Toc() / 100; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Tic();
  cdouble UD2;
  for(int n=0; n<100; n++)
   UD2 = GetGMatrixElement_DA(O1, nfa, O2, nfb, Omega, 2); 
  double TimeD2 = Toc() / 100; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble URef = RefSI ? USI : UVI;
   
  printf("--nfa %i --nfb %i\n",nfa,nfb);
  printf("URef=(%e,%e)\n",real(URef),imag(URef));
  printf("VI(16)  RD %.1e, Time %e ms \n",RD(UVI16, URef), TimeVI16*1.0e3);
  printf("SI(0):  RD %.1e, Time %e ms \n",RD(USI,   URef),  Time0*1.0e3);
  printf("SI(5):  RD %.1e, Time %e ms \n",RD(USI5,  URef),  Time5*1.0e3);
  printf("SI(9):  RD %.1e, Time %e ms \n",RD(USI9,  URef),  Time9*1.0e3);
  printf("SI(16): RD %.1e, Time %e ms \n",RD(USI16, URef), Time16*1.0e3);
  printf("SI(20): RD %.1e, Time %e ms \n",RD(USI20, URef), Time20*1.0e3);
  printf("D1:     RD %.1e, Time %e ms \n",RD(UD1,   URef),   TimeD1*1.0e3);
  printf("D2:     RD %.1e, Time %e ms \n",RD(UD2,   URef),   TimeD2*1.0e3);

}
