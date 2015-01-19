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

namespace buff { 

int GetVInverseElements(SWGVolume *V, int nfA, cdouble Omega,
                        int Indices[7], cdouble Entries[7]);

int GetOverlapElements(SWGVolume *O, int nfA,
                       int Indices[7], double Entries[7]);

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
  cdouble Omega=1.0;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,          0, "frequency"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry(GeoFile);
  SWGVolume *O   = G->Objects[0];
 
  srand48(time(0));
  if (nfA==-1)
   nfA = lrand48() % O->NumInteriorFaces;
  
  SWGFace *F = O->Faces[nfA];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int Indices1[7]; 
  double Entries1[7];
  int NNZ1  = GetOverlapElements(O, nfA, Indices1, Entries1);

  int Indices2[7];
  cdouble Entries2[7];
  int NNZ2 = GetVInverseElements(O, nfA, Omega, Indices2, Entries2);

  if (NNZ1!=NNZ2)
   ErrExit("bawonkatage! NNZ1=%i, NNZ2=%i!\n",NNZ1,NNZ2);

  cdouble Eps[3][3];
  O->MP->GetEps(Omega, F->Centroid, Eps);
  double PreFac = real( 1.0 / (Omega*Omega*(1.0-Eps[0][0])) );

  for(int nnz=0; nnz<NNZ1; nnz++)
   printf("<%3i|%3i (%3i)>: (%+.3e) (%+.3e,%+.1e) (%.1e) (%.2e)\n",
           nfA,Indices1[nnz],Indices2[nnz],
           PreFac*Entries1[nnz],
           real(Entries2[nnz]), imag(Entries2[nnz]),
           RD(PreFac*Entries1[nnz], Entries2[nnz]),
           PreFac*Entries1[nnz] / real(Entries2[nnz])
         );

}
