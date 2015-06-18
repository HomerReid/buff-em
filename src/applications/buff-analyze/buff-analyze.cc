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
/* quality factor for tetrahedron, defined as                  */
/*  volume / ( (avg edge length) * (total surface area) )      */
/* normalized by the quality factor for an ideal tetrahedron   */
/***************************************************************/
double GetTetrahedronQualityFactor(SWGVolume *O, int nt)
{
  SWGTet *T = O->Tets[nt];
  
  double TotalArea =  O->Faces[T->FI[0]]->Area
                     +O->Faces[T->FI[1]]->Area
                     +O->Faces[T->FI[2]]->Area
                     +O->Faces[T->FI[3]]->Area;

  double *V1 = O->Vertices + 3*T->VI[0];
  double *V2 = O->Vertices + 3*T->VI[1];
  double *V3 = O->Vertices + 3*T->VI[2];
  double *V4 = O->Vertices + 3*T->VI[3];
  double AvgLength =
   ( VecDistance(V2,V1) + VecDistance(V3,V1) + VecDistance(V4,V1)
    +VecDistance(V3,V2) + VecDistance(V4,V2) + VecDistance(V4,V3)
   ) / 6.0;

  // quality factor for equilateral tetrahedron
  #define IDEAL_QUALITY_FACTOR 1.0/(6.0*sqrt(6.0))

  return T->Volume / (TotalArea*AvgLength*IDEAL_QUALITY_FACTOR);

}

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

  printf(" Volume:        %e \n",Volume);
  printf(" Surface area:  %e \n",SurfaceArea);
  printf(" \n");

  double MinQF=1.0e9, MaxQF=-1.0e9, AvgQF=0.0;
//  FILE *f=vfopen("/tmp/%s.TQFs","w",GetFileBase(V->MeshFileName));
  for(int nt=0; nt<V->NumTets; nt++)
   { double QF = GetTetrahedronQualityFactor(V,nt);
#if 0
     fprintf(f,"Tet %5i: (V,QF) = (%e,%e) A=(%.2e,%.2e,%.2e,%.2e)\n",
                nt,V->Tets[nt]->Volume,QF,
                V->Faces[V->Tets[nt]->FI[0]]->Area,
                V->Faces[V->Tets[nt]->FI[1]]->Area,
                V->Faces[V->Tets[nt]->FI[2]]->Area,
                V->Faces[V->Tets[nt]->FI[3]]->Area);
#endif
     MinQF=fmin(MinQF, QF);
     MaxQF=fmax(MaxQF, QF);
     AvgQF+=QF;
   };
//  fclose(f);
  AvgQF/=((double)(V->NumTets));

  printf(" Min/Max/Avg Quality factor: %e / %e / %e \n",
          MinQF, MaxQF, AvgQF);
  printf("\n");
}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* convenience shortcuts that allow the code to be invoked as  */
  /*  % buff-analyze File.vmsh                                   */
  /* or                                                          */
  /*  % buff-analyze File.buffgeo                                */
  /***************************************************************/
  char *GeoFile=0;
  char *MeshFile=0;
  if (argc>=2)
   { char *Ext = GetFileExtension(argv[1]);
     if (Ext && !StrCaseCmp(Ext,"buffgeo"))
      { GeoFile=strdup(argv[1]); 
        argv[1]=0; 
      }
     else if (Ext && (!StrCaseCmp(Ext,"msh") || !StrCaseCmp(Ext,"vmsh") ) )
      { MeshFile=strdup(argv[1]);
        argv[1]=0; 
      };
   };

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *TransFile=0;
  bool WriteGMSHFiles=false;
  bool PlotPermittivity=false;
  int PlotBF[10], nPlotBF;
  int PlotTet[10], nPlotTet;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING, 1, 1, (void *)&GeoFile,          0, ".buffgeo file"},
     {"mesh",               PA_STRING, 1, 1, (void *)&MeshFile,         0, ".msh file"},
     {"meshfile",           PA_STRING, 1, 1, (void *)&MeshFile,         0, ".msh file"},
     {"transfile",          PA_STRING, 1, 1, (void *)&TransFile,        0, "list of geometrical transformations"},
     {"WriteGMSHFiles",     PA_BOOL,   0, 1, (void *)&WriteGMSHFiles,   0, "output GMSH visualization code"},
     {"PlotPermittivity",   PA_BOOL,   0, 1, (void *)&PlotPermittivity, 0, "color tetrahedra by average local permittivity"},
     {"PlotBF",             PA_INT,    1, 10, (void *)PlotBF,    &nPlotBF, "draw an individual basis function"},
     {"PlotTet",            PA_INT,    1, 10, (void *)PlotTet,   &nPlotTet, "draw an individual tetrahedron"},
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

  /***************************************************************/
  /* plot individual basis functions if that was requested       */
  /***************************************************************/
  if (nPlotBF)
   { FILE *f=vfopen("%s.BFs.pp","w",GetFileBase(GeoFile));
     for(int n=0; n<nPlotBF; n++)
      G->Objects[0]->DrawBF(PlotBF[n], f);
     fclose(f);
     printf("Individual BF visualizations written to file %s.BFs.pp",
             GetFileBase(GeoFile));
   };

  /***************************************************************/
  /* plot individual tetrahedra if that was requested            */
  /***************************************************************/
  if (nPlotTet)
   { FILE *f=vfopen("%s.Tets.pp","w",GetFileBase(GeoFile));
     for(int n=0; n<nPlotTet; n++)
      G->Objects[0]->DrawTet(PlotTet[n], f);
     fclose(f);
     printf("Individual tetrahedron visualizations written to file %s.BFs.pp",
             GetFileBase(GeoFile));
   };

  /***************************************************************/
  /* plot the overall geometry if that was requested             */
  /***************************************************************/
  if (WriteGMSHFiles && GeoFile!=0)
   { char FileName[1000];
     sprintf(FileName,"%s.pp",GetFileBase(GeoFile));
     G->WritePPMesh(FileName,"Default");
     printf("GMSH visualization data written to %s.\n",FileName);
   };

  /***************************************************************/
  /* plot the permittivity if that was requested                 */
  /***************************************************************/
  if (PlotPermittivity && GeoFile!=0)
   { char FileName[1000];
     sprintf(FileName,"%s.Epsilon.pp",GetFileBase(GeoFile));
     G->PlotPermittivity(FileName,"Default");
     printf("GMSH permittivity plot written to %s.\n",FileName);
   };

  /***************************************************************/
  /* handle geometrical transformations if any were specified    */
  /***************************************************************/
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
