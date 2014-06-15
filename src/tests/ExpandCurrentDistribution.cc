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
  int nfA=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,        0, "mesh file"},
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
  /* form the overlap matrix *************************************/
  /***************************************************************/
  HMatrix *M = G->AllocateVIEMatrix();
  M->Zero();
  for(int nfa=0; nfa<O->NumInteriorFaces; nfa++)
   { int Indices[7];
     cdouble Entries[7];
     int NNZ = GetVInverseElements(O, nf, Omega, Indices, Entries);
     for(int nnz=0; nnz<NNZ; nnz++)
      { int nfb = Indices[nnz];
        M->SetEntry(RowOffset + nfa, ColOffset + nfb, Entries[nnz]);
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *J = G->AllocateRHSVector();

}
