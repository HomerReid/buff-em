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
 * tNSIntegrals.cc
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

cdouble GetGMatrixElement_VI(SWGVolume *VA, int nfA,
                             SWGVolume *VB, int nfB,
                             cdouble Omega, int Order);
              }

#define II cdouble(0.0,1.0)

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
  int nfA=-1;
  bool SingularCase=false;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"nfA",                PA_INT,     1, 1, (void *)&nfA,            0, "nfA"},
     {"SingularCase",       PA_BOOL,    0, 1, (void *)&SingularCase,   0, "SingularCase"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry(GeoFile);

  SWGVolume *OA = G->Objects[0];

  srand48(time(0));
  if (nfA==-1)
   nfA = lrand48() % OA->NumInteriorFaces;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen("tAllIntegrals.preamble","w");
  fprintf(f,"#columns: \n");
  fprintf(f,"#1,2   = nfA, nfB\n");
  fprintf(f,"#3     = number of common vertices \n");
  fprintf(f,"#4     = relative distance \n");
  fprintf(f,"#5 6   = real, imag U\n"); 
  fprintf(f,"#7,8   = RD, time for VI4\n");
  fprintf(f,"#9,10  = RD, time for VI16\n");
  fprintf(f,"#11,12 = RD, time for SI5\n");
  fprintf(f,"#13,14 = RD, time for SI9\n");
  fprintf(f,"#15,16 = RD, time for SI16\n");
  fprintf(f,"#17,18 = RD, time for D1\n");
  fprintf(f,"#19,20 = RD, time for D2\n");
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  f=fopen("tAllIntegrals.out","w");
  for(int no=0; no<G->NumObjects; no++)
   for(int nfB=0; nfB<G->Objects[no]->NumInteriorFaces; nfB++)
    { 
      SWGVolume *OB = G->Objects[no];
      double rRel;
      int ncv = CompareBFs(OA, nfA, OB, nfB, &rRel);

      cdouble URef;
      if (SingularCase)
       { if (ncv==0) continue;
         URef=GetGMatrixElement_SI(OA, nfA, OB, nfB, Omega, 0);
       }
      else
       { if (ncv!=0) continue;
         URef=GetGMatrixElement_VI(OA, nfA, OB, nfB, Omega, 0);
       };

      Tic();
      cdouble UVI4=GetGMatrixElement_VI(OA, nfA, OB, nfB, Omega, 4);
      double TVI4=Toc();

      Tic();
      cdouble UVI16=GetGMatrixElement_VI(OA, nfA, OB, nfB, Omega, 16);
      double TVI16=Toc();

      Tic();
      cdouble USI5 = GetGMatrixElement_SI(OA, nfA, OB, nfB, Omega, 5);
      double TSI5 = Toc();

      Tic();
      cdouble USI9 = GetGMatrixElement_SI(OA, nfA, OB, nfB, Omega, 9);
      double TSI9 = Toc();

      Tic();
      cdouble USI16 = GetGMatrixElement_SI(OA, nfA, OB, nfB, Omega, 16);
      double TSI16 = Toc();

      Tic();
      cdouble USI20 = GetGMatrixElement_SI(OA, nfA, OB, nfB, Omega, 20);
      double TSI20 = Toc();

      Tic();
      cdouble UD1 = GetGMatrixElement_DA(OA, nfA, OB, nfB, Omega, 1); 
      double TD1 = Toc();

      Tic();
      cdouble UD2 = GetGMatrixElement_DA(OA, nfA, OB, nfB, Omega, 2); 
      double TD2 = Toc();

      printf("%i %i %i %e\n",nfA,nfB,ncv,rRel);
      fprintf(f,"%i %i %i %e ",nfA,nfB,ncv,rRel); 
      fprintf(f,"%e %e ",real(URef),imag(URef));
      fprintf(f,"%e %e ",RD(URef,UVI4),  TVI4);
      fprintf(f,"%e %e ",RD(URef,UVI16), TVI16);
      fprintf(f,"%e %e ",RD(URef,USI5),  TSI5);
      fprintf(f,"%e %e ",RD(URef,USI9),  TSI9);
      fprintf(f,"%e %e ",RD(URef,USI16), TSI16);
      fprintf(f,"%e %e ",RD(URef,USI20), TSI20);
      fprintf(f,"%e %e ",RD(URef,UD1),   TD1);
      fprintf(f,"%e %e ",RD(URef,UD2),   TD2);
      fprintf(f,"\n");

    };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int Status;
  Status=system("sort -g -k4 tAllIntegrals.out > tAllIntegrals.out.sorted");
  Status=system("sort -g -k4 tAllIntegrals.out > tAllIntegrals.out.sorted");
  Status=system("cat tAllIntegrals.preamble tAllIntegrals.out.sorted > tAllIntegrals.dat");

}
