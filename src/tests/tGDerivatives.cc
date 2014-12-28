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
 * homer reid       -- 8/2014
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
                             cdouble Omega, 
                             int Order=0, cdouble *dG=0);

cdouble GetGMatrixElement_VI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, int Order,
                             cdouble *GradG);
              }

namespace scuff{
void CalcGC(double R[3], cdouble Omega,
            cdouble EpsR, cdouble MuR,
            cdouble GMuNu[3][3], cdouble CMuNu[3][3],
            cdouble GMuNuRho[3][3][3], cdouble CMuNuRho[3][3][3]);
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

  cdouble k  = Data->k;
  double Eta = Data->Eta;

  double R[3];
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);

  double XmX0[3];
  XmX0[0] = xA[0] - 0.0;
  XmX0[1] = xA[1] - 0.0;
  XmX0[2] = xA[2] - 0.0;

  cdouble GMuNu[3][3], GMuNuRho[3][3][3];
  CalcGC(R, k, 1.0, 1.0, GMuNu, 0, GMuNuRho, 0);

  cdouble GMuNuTheta[3][3][3];
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    for(int Rho=0; Rho<3; Rho++)
     {
       int RP1=(Rho+1)%3, RP2=(Rho+2)%3;
       GMuNuTheta[Mu][Nu][Rho]
        = (XmX0[RP1]*GMuNuTheta[Mu][Nu][RP2] - XmX0[RP2]*GMuNuTheta[Mu][Nu][RP1]);
     };

  cdouble *zF=(cdouble *)I;
  memset(zF, 0, 7*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { zF[0] += bA[Mu] * GMuNu[Mu][Nu]         * bB[Nu];
      zF[1] += bA[Mu] * GMuNuRho[Mu][Nu][0]   * bB[Nu];
      zF[2] += bA[Mu] * GMuNuRho[Mu][Nu][1]   * bB[Nu];
      zF[3] += bA[Mu] * GMuNuRho[Mu][Nu][2]   * bB[Nu];
      zF[4] += bA[Mu] * GMuNuTheta[Mu][Nu][0] * bB[Nu];
      zF[5] += bA[Mu] * GMuNuTheta[Mu][Nu][1] * bB[Nu];
      zF[6] += bA[Mu] * GMuNuTheta[Mu][Nu][2] * bB[Nu];
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFDDerivatives(SWGVolume *OACopy, int nfa, SWGVolume *OB, int nfb,
                      cdouble Omega, cdouble G0, cdouble dG[6])
{
  
  for(int Mu=0; Mu<3; Mu++)
   { 
     double DX[3] = {0.0, 0.0, 0.0};

     DX[Mu] = 0.001;
     OACopy->Transform("DISPLACED %e %e %e \n",DX[0],DX[1],DX[2]);
     cdouble GpdG = GetGMatrixElement_SI(OACopy, nfa, OB, nfb, Omega, 20);
     OACopy->UnTransform();

     DX[Mu] = -0.001;
     OACopy->Transform("DISPLACED %e %e %e \n",DX[0],DX[1],DX[2]);
     cdouble GmdG = GetGMatrixElement_SI(OACopy, nfa, OB, nfb, Omega, 20);
     OACopy->UnTransform();

     dG[Mu] = (GpdG - GmdG) / (-2.0*DX[Mu]);
   };

  for(int Mu=0; Mu<3; Mu++)
   { 
     double nHat[3] = {0.0, 0.0, 0.0};
     nHat[Mu] = 1.0;

     double DeltaTheta = 0.01; // radians 
     OACopy->Transform("ROTATED %e ABOUT %e %e %e ",DeltaTheta*180.0/M_PI,nHat[0],nHat[1],nHat[2]);
     cdouble GpdG = GetGMatrixElement_SI(OACopy, nfa, OB, nfb, Omega, 20);
     OACopy->UnTransform();

     DeltaTheta = -0.01;
     OACopy->Transform("ROTATED %e ABOUT %e %e %e ",DeltaTheta*180.0/M_PI,nHat[0],nHat[1],nHat[2]);
     cdouble GmdG = GetGMatrixElement_SI(OACopy, nfa, OB, nfb, Omega, 20);
     OACopy->UnTransform();

     dG[3+Mu] = (GpdG - GmdG) / (-2.0*DeltaTheta);
   };

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
  bool DoBF=false;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nfa",                PA_INT,     1, 1, (void *)&nfa,            0, "nfa"},
     {"nfb",                PA_INT,     1, 1, (void *)&nfb,            0, "nfb"},
     {"ncv",                PA_INT,     1, 1, (void *)&ncv,            0, "ncv"},
     {"DoBF",               PA_BOOL,    0, 1, (void *)&DoBF,           0, "do the brute-force calculation"},
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
  cdouble GBF[7];
  if (DoBF)
   { cdouble Error[7];
     GVIData MyData, *Data=&MyData;
     Data->k=Omega;
     Data->Eta = 0.0;
     BFBFInt(O, nfa, O, nfb, GVIntegrand, (void *)Data,
             14, (double *)GBF, (double *)Error, 0, 1000000, 1.0e-6);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble GVI[7];
  GVI[0] = GetGMatrixElement_VI(O, nfa, O, nfb, Omega, 16, GVI+1);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble GSI[7];
  memset(GSI, 0, 7*sizeof(cdouble));
  GSI[0] = GetGMatrixElement_SI(O, nfa, O, nfb, Omega, 20, GSI+1);

  /***************************************************************/
  /* If we didn't do the BF calculation then estimate derivatives*/
  /* by finite-differencing                                      */
  /***************************************************************/
  if (!DoBF)
   { 
     SWGGeometry *GCopy = new SWGGeometry(GeoFile);
     SWGVolume *OCopy = GCopy->Objects[0];
     GetFDDerivatives(OCopy, nfa, O, nfb, Omega, GSI[0], GBF+1);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetDefaultCD2SFormat("(%+.2e, %+.2e)");
  printf(" n | %s | %s |   RD    | %s |   RD    \n",
         "         BF           ",
         "         VI           ",
         "         SI           ");
  for(int nf=0; nf<7; nf++)
   { printf(" %i | %s | %s | %.1e | %s |%.1e\n", nf, CD2S(GBF[nf]), 
            CD2S(GVI[nf]), RD(GBF[nf],GVI[nf]),
            CD2S(GSI[nf]), RD(GBF[nf],GSI[nf]));
   };

}
