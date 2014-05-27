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
 * TGFluctuations.cc -- a code for predicting fluctuation phenomena
 *                   -- (Casimir forces and radiative heat transfer) in
 *                   -- 3D bodies with inhomogeneous anisotropic 
 *                   -- material properties
 *
 * homer reid       -- 5/2014
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
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *MeshFile1=0;
  char *MatFile1=0;
  double T1=0.0;
  char *MeshFile2=0;
  char *MatFile2=0;
  double T2=0.0;
  char *TransFile=0;
  char *OmegaFile=0;

  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"mesh1",       PA_STRING,  1, 1, (void *)&MeshFile1, 0, "mesh file 1"},
     {"material1",   PA_STRING,  1, 1, (void *)&MatFile1,  0, "material file 1"},
     {"temp1"        PA_STRING,  1, 1, (void *)&T1,        0, "temperature 1 in Kelvin"},
//
     {"mesh2",       PA_STRING,  1, 1, (void *)&MeshFile2, 0, "mesh file 2"},
     {"material2",   PA_STRING,  1, 1, (void *)&MatFile2,  0, "material file 2"},
     {"temp2"        PA_STRING,  1, 1, (void *)&T2,        0, "temperature 2 in Kelvin"},
//
     {"TransFile",   PA_STRING,  1, 1, (void *)&TransFile, 0, "list of geometrical transformations"},
//
     {"OmegaFile",   PA_STRING,  1, 1, (void *)&OmegaFile, 0, "list of Omega values"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (MeshFile1==0 || MatFile1)
   OSUsage(argv[0],OSArray,"--mesh1 and --material1 options are mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGVolume *V1    = new SWGVolume(MeshFile1);
  IHAIMatProp *MP1 = new IHAIMatProp(MatFile1);

  SWGVolume *V2    = MeshFile2 ? new SWGVolume(MeshFile2)  : 0;
  IHAIMatProp *MP2 = MatFile2  ? new IHAIMatProp(MatFile2) : 0;
  if ( V2!=0 && M2==0 || V2==0 && M2!=0 )
   OSUsage(argv[0],OSArray,"--mesh2 and --material2 options must be both present or both absent");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HVector *OmegaList = OmegaFile ? new HVector(OmegaFile) : 0;
  if (OmegaList && OmegaList->ErrMsg)
   ErrExit(OmegaList->ErrMsg);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if Integrands
  int NumIntegrands = (V2) ?

  for(int nOmega=0; nOmega<OmegaList->N; nOmega++)
   { GetIntegrands(V1, MP1, V2, MP2, Omega, Integrands);
     GetIntegrands(V1, MP1, V2, MP2, Omega);

}
