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
 * homer reid -- 6/2014
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

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble Omega=0.1;
  char *EPFile=0;
  int nfA=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "omega"},
     {"EPFile",             PA_STRING,  1, 1, (void *)&EPFile,         0, "mesh file"},
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  // E = (Eps-1)/(Eps+2) * EInc
  // J = (Eps-1)/(Eps+2) * EInc
  cdouble Eps[3][3], Eps00;
  double X[3]={0.0, 0.0, 0.0};
  G->Objects[0]->MP->GetEps(Omega, X, Eps);
  Eps00=Eps[0][0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *J = G->AllocateRHSVector();
  double pwDir[3]={0.0, 0.0, 1.0};
  cdouble pwPol[3]={1.0, 0.0, 0.0};
  pwPol[0] = (3.0 / (Eps00 + 2.0)) * (-II*Omega*ZVAC*(Eps00-1.0));
  PlaneWave *PW = new PlaneWave(pwPol, pwDir);
  G->AssembleRHSVector(0.01, PW, J);
  G->PlotCurrentDistribution("ECD.pp",J,"RHS");
 
  /***************************************************************/
  /* form the overlap matrix *************************************/
  /***************************************************************/
  HMatrix *M = G->AllocateVIEMatrix();
  M->Zero();
  for(int nfa=0; nfa<O->NumInteriorFaces; nfa++)
   { int Indices[7];
     double Entries[7];
     int NNZ = GetOverlapElements(O, nfa, Indices, Entries);
     for(int nnz=0; nnz<NNZ; nnz++)
      { int nfb = Indices[nnz];
        M->SetEntry(nfa, nfb, Entries[nnz]);
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  M->LUFactorize();
  M->LUSolve(J);
  G->PlotCurrentDistribution("ECD.pp",J,"J");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *XMatrix=new HMatrix(EPFile,LHM_TEXT);
  if (XMatrix->ErrMsg)
   ErrExit("Error processing EP file: %s\n",XMatrix->ErrMsg);

  HMatrix *FMatrix1 = G->GetFields( 0, J, Omega, XMatrix); // scattered
  HMatrix *FMatrix2 = G->GetFields(PW, 0, Omega, XMatrix); // incident

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char buffer[MAXSTR];
  snprintf(buffer,MAXSTR,"%s.scattered",GetFileBase(EPFile));
  FILE *f1=CreateUniqueFile(buffer,1);
  snprintf(buffer,MAXSTR,"%s.total",GetFileBase(EPFile));
  FILE *f2=CreateUniqueFile(buffer,1);

  int nr, nc; 
  SetDefaultCD2SFormat("%.8e %.8e ");
  for(nr=0; nr<FMatrix1->NR; nr++)
   { fprintf(f1,"%.8e %.8e %.8e ",XMatrix->GetEntryD(nr, 0),
                                  XMatrix->GetEntryD(nr, 1),
                                  XMatrix->GetEntryD(nr, 2));

     fprintf(f2,"%.8e %.8e %.8e ",XMatrix->GetEntryD(nr, 0),
                                  XMatrix->GetEntryD(nr, 1),
                                  XMatrix->GetEntryD(nr, 2));

     for(nc=0; nc<FMatrix1->NC; nc++)
      { 
        fprintf(f1,"%s ",CD2S(  FMatrix1->GetEntry(nr,nc)) );

        fprintf(f2,"%s ",CD2S(  FMatrix1->GetEntry(nr,nc)  
                               +FMatrix2->GetEntry(nr,nc)) );
      };

     fprintf(f1,"\n");
     fprintf(f2,"\n");

   };
  fclose(f1);
  fclose(f2);

}
