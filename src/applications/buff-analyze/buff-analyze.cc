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

#define FIBBIDATALEN 6

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
  /* MOI=moment of inertia ***************************************/
  /***************************************************************/
  double Volume=0.0, SurfaceArea=0.0, MOI[3]={0.0,0.0,0.0};
  for(int nt=0; nt<V->NumTets; nt++)
   { SWGTet *T=V->Tets[nt];
     double *X0=V->Tets[nt]->Centroid;
     Volume+=T->Volume;
     MOI[0]+=T->Volume*(X0[1]*X0[1] + X0[2]*X0[2]);
     MOI[1]+=T->Volume*(X0[2]*X0[2] + X0[0]*X0[0]);
     MOI[2]+=T->Volume*(X0[0]*X0[0] + X0[1]*X0[1]);
   };
  for(int nf=V->NumInteriorFaces; nf<V->NumTotalFaces; nf++)
   SurfaceArea+=V->Faces[nf]->Area;

  printf(" Volume:          %e \n",Volume);
  printf(" Surface area:    %e \n",SurfaceArea);
  printf(" \n");
  printf(" Moment of inertia: {%+.2e,%+.2e,%+.2e}\n",MOI[0],MOI[1],MOI[2]);
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
/***************************************************************/
/***************************************************************/
void WriteCache(SWGVolume *O, int NumChunks=1, int WhichChunk=0)
{
  FIBBICache *Cache = new FIBBICache(O->MeshFileName);

  int NF     = O->NumInteriorFaces;
  Log("Writing GCache for %s", O->MeshFileName);
  if (NumChunks>1)
   LogC(" (chunk %i/%i)",WhichChunk,NumChunks);

  unsigned long M0=GetMemoryUsage() / (1<<20);
  Log("Initial memory usage: %8lu MB",M0);

  int NumPairs = NF*(NF+1)/2;
  int ChunkSize = NumPairs / NumChunks;
  int npMin = WhichChunk*ChunkSize;
  int npMax = npMin + ChunkSize - 1;
  if (npMax >= NumPairs)
   npMax=NumPairs-1;
  int np=-1;
  int NumRecords=0;
#ifdef USE_OPENMP
  int NumThreads=GetNumThreads();
  Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      \
                         reduction(+:NumRecords)   \
                         num_threads(NumThreads)
#endif
  for(int nfa=0; nfa<NF; nfa++)
   for(int nfb=nfa; nfb<NF; nfb++)
    {
      //if (nfb<nfa) continue;

      np++;
      if (np<npMin) continue;
      if (np>=npMax) continue;

      int ncv = CompareBFs(O, nfa, O, nfb);
      if (ncv==0) continue;

      LogPercent(np-npMin, ChunkSize, 100);
      double Data[FIBBIDATALEN];
      Cache->GetFIBBIData(O, nfa, O, nfb, Data);
      NumRecords++;
    };
// end of multithreaded loop

  char FileName[MAXSTR];
  if (NumChunks>1)
   {
     snprintf(FileName,MAXSTR,"%s.cache.%i.%i",
                               RemoveExtension(O->MeshFileName),
                               WhichChunk, NumChunks);
   }
  else
   {
     snprintf(FileName,MAXSTR,"%s.cache",RemoveExtension(O->MeshFileName));
   };
  Cache->Store(FileName);

  printf("Wrote %i FIBBI records to %s.\n",NumRecords,FileName);

  unsigned long M1=GetMemoryUsage() / (1<<20);
  Log("Memory with cache:    %8lu MB (cache size %8lu) ",M1,M1-M0);
  delete Cache;

  unsigned long M2=GetMemoryUsage() / (1<<20);
  Log("Memory after free:    %8lu MB (leaked:    %8lu)",M2,M2-M0);

}

/***************************************************************/
/* report permittivity at user-specified (Omega, X, Y, Z)      */
/***************************************************************/
void GetPermittivity(char *GeoFile, cdouble Omega, double XYZ[3])
{
   SWGGeometry *G = new SWGGeometry(GeoFile);
   printf("Permittivities at (Omega, x, y, z)=(%g,%g,%g,%g): \n",
           real(Omega),XYZ[0],XYZ[1],XYZ[2]);

   FILE *f=vfopen("%s.permittivity","r",GetFileBase(GeoFile));
   if (f==0)
    { f=vfopen("%s.permittivity","w",GetFileBase(GeoFile));
      fprintf(f,"# data file columns: \n");
      fprintf(f,"# 1     omega \n");
      fprintf(f,"# 2,3,4 x,y,z \n");
      fprintf(f,"# 5     index of object in .buffgeo file\n");
      fprintf(f,"# 6,7   re,im EpsXX\n");
      fprintf(f,"# 8,9   re,im EpsXY\n");
      fprintf(f,"# 10,11 re,im EpsXZ\n");
      fprintf(f,"# 12,13 re,im EpsYX\n");
      fprintf(f,"# 14,15 re,im EpsYY\n");
      fprintf(f,"# 16,17 re,im EpsYZ\n");
      fprintf(f,"# 18,19 re,im EpsZX\n");
      fprintf(f,"# 20,21 re,im EpsZY\n");
      fprintf(f,"# 22,23 re,im EpsZZ\n");
      fclose(f);
    };

   f=vfopen("%s.permittivity","a",GetFileBase(GeoFile));
   for(int no=0; no<G->NumObjects; no++)
    { 
      SWGVolume *O=G->Objects[no];
      SVTensor *SVT=O->SVT;
      cdouble Eps[3][3];
      SVT->Evaluate(Omega, XYZ, Eps);

      SetDefaultCD2SFormat("{%+.3e, %+.3e}");
      printf("For object %s (%s): ",O->Label,O->MeshFileName);
      if (SVT->Isotropic)
       printf(" %s \n",CD2S(Eps[0][0]));
      else
       { printf("\n");
         printf(" %s %s %s \n",CD2S(Eps[0][0]),CD2S(Eps[0][1]),CD2S(Eps[0][2]));
         printf(" %s %s %s \n",CD2S(Eps[1][0]),CD2S(Eps[1][1]),CD2S(Eps[1][2]));
         printf(" %s %s %s \n",CD2S(Eps[2][0]),CD2S(Eps[2][1]),CD2S(Eps[2][2]));
         printf("\n");
       };
      printf("\n");

      SetDefaultCD2SFormat("%e %e");
      fprintf(f,"%e %e %e %e %i ",real(Omega),XYZ[0],XYZ[1],XYZ[2],no);
      for(int Mu=0; Mu<3; Mu++)
       for(int Nu=0; Nu<3; Nu++)
        fprintf(f,"%s ",CD2S(Eps[Mu][Nu]));
      fprintf(f,"\n");
    };

   fclose(f);
   printf("Thank you for your support.\n");
   exit(0);

}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  SetLogFileName("buff-analyze.log");

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
  bool WriteGCache=false;
  int NumChunks=1;
  int WhichChunk=0;
  double XYZ[3]; int nXYZ=0;
  cdouble Omega=1.0;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1,  (void *)&GeoFile,          0, ".buffgeo file"},
     {"mesh",               PA_STRING,  1, 1,  (void *)&MeshFile,         0, ".msh file"},
     {"meshfile",           PA_STRING,  1, 1,  (void *)&MeshFile,         0, ".msh file"},
/**/
     {"transfile",          PA_STRING,  1, 1,  (void *)&TransFile,        0, "list of geometrical transformations"},
/**/
     {"WriteGMSHFiles",     PA_BOOL,    0, 1,  (void *)&WriteGMSHFiles,   0, "output GMSH visualization code"},
/**/
     {"XYZ",                PA_DOUBLE,  3, 1,  (void *)XYZ,          &nXYZ,  "coordinates at which to evaluate permittivity"},
     {"Omega",              PA_CDOUBLE, 1, 1, (void *)&Omega,             0, "angular frequency for permittivity check"},
/**/
     {"WriteCache",         PA_BOOL,    0, 1, (void *)&WriteGCache,       0, "write cache file"},
     {"WriteGCache",        PA_BOOL,    0, 1, (void *)&WriteGCache,       0, "write cache file"},
     {"NumChunks",          PA_INT,     1, 1, (void *)&NumChunks,         0, "number of pieces into which to subdivide cache write (1)"},
     {"WhichChunk",         PA_INT,     1, 1, (void *)&WhichChunk,        0, "which piece to write (0)"},
/**/
     {"PlotPermittivity",   PA_BOOL,    0, 1, (void *)&PlotPermittivity,  0, "color tetrahedra by average local permittivity"},
     {"PlotBF",             PA_INT,     1, 10, (void *)PlotBF,     &nPlotBF, "draw an individual basis function"},
     {"PlotTet",            PA_INT,     1, 10, (void *)PlotTet,   &nPlotTet, "draw an individual tetrahedron"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (GeoFile && nXYZ)
   GetPermittivity(GeoFile, Omega, XYZ);

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
   { SWGVolume *O=new SWGVolume(MeshFile);
     if (WriteGCache)
      WriteCache(O, NumChunks, WhichChunk);
     else
      AnalyzeVolume( O );
   }
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

  printf("Thank you for your support.\n");

}
