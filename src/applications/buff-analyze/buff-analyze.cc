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
 * buff-analyze.cc
 *
 * homer reid       -- 4/2014
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
void AnalyzeVolume(SWGVolume *V)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Volume %s (file %s): \n",V->Label,V->MeshFileName);
  printf(" %3i vertices \n",V->NumVertices);
  printf(" %3i tetrahedra \n",V->NumTets);
  printf(" %3i interior faces \n",V->NumInteriorFaces);
  printf(" %3i exterior faces \n",V->NumTotalFaces - V->NumInteriorFaces);
  printf("\n");
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Volume=0.0, SurfaceArea=0.0;
  for(int nt=0; nt<V->NumTets; nt++)
   Volume+=V->Tets[nt]->Volume;
  for(int nf=V->NumInteriorFaces; nf<V->NumTotalFaces; nf++)
   SurfaceArea+=V->Faces[nf]->Area;

  printf(" Volume:       %e \n",Volume);
  printf(" Surface area: %e \n",SurfaceArea);
  printf("\n");
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
  char *MeshFile=0;
  char *TransFile=0;
  bool WriteGMSHFiles=false;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING, 1, 1, (void *)&GeoFile,        0, "mesh file"},
     {"mesh",               PA_STRING, 1, 1, (void *)&MeshFile,       0, "mesh file"},
     {"meshfile",           PA_STRING, 1, 1, (void *)&MeshFile,       0, "mesh file"},
     {"WriteGMSHFiles",     PA_BOOL,   0, 1, (void *)&WriteGMSHFiles, 0, "output GMSH visualization code"},
     {"transfile",          PA_STRING, 1, 1, (void *)&TransFile,      0, "output GMSH visualization code"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G=0;
  if (GeoFile)
   { 
     G = new SWGGeometry(GeoFile);
     printf("\n*\n* Geometry file %s: %i objects \n*\n\n",
             GeoFile,G->NumObjects);
     for(int no=0; no<G->NumObjects; no++)
      { printf("**************************************************\n");
        printf("* Object %i: \n",no);
        printf("**************************************************\n");
        if ( G->Mate[no] != -1 )
         printf(" (duplicate of object %i, %s)\n\n",
                  G->Mate[no],
                  G->Objects[G->Mate[no]]->Label);
        else
         AnalyzeVolume( G->Objects[no] );
      };
   }
  else if (MeshFile)
   AnalyzeVolume( new SWGVolume(MeshFile) );
  else
   ErrExit("either --geometry or --mesh / --meshfile option is mandatory");

  if (WriteGMSHFiles && GeoFile!=0)
   { char FileName[1000];
     sprintf(FileName,"%s.pp",GetFileBase(GeoFile));
     G->WritePPMesh(FileName,"Default"); 
     printf("GMSH visualization data written to %s.\n",FileName);
   };

  if (TransFile)
   { 
     if (!G) 
      ErrExit("--transfile may only be used with --geometry");

     int ngtc, NGTC;
     GTComplex **GTCList=ReadTransFile(TransFile, &NGTC);
     char *ErrMsg=G->CheckGTCList(GTCList, NGTC);
     if (ErrMsg)
      ErrExit("file %s: %s",TransFile,ErrMsg);

     char FileName[1000];
     sprintf(FileName,"%s.transformed.pp",GetFileBase(G->GeoFileName));
     unlink(FileName);

     for(ngtc=0; ngtc<NGTC; ngtc++) 
      {
        G->Transform(GTCList[ngtc]);
        G->WritePPMesh(FileName, GTCList[ngtc]->Tag);
        G->UnTransform();
      };

     printf("Visualizations for %i transforms written to %s.\n",
             NGTC,FileName);
   };

}
