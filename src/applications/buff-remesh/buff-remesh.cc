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
 * buff-remesh.cc
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
double GetTetQF(double *V[4], double *VV, double *AA, double *LL, FILE *f)
{ 
  double A[3], B[3], C[3], D[3], E[3], F[3], Scratch[3];

  VecSub(V[1],V[0],A);
  VecSub(V[2],V[0],B);
  VecSub(V[3],V[0],C);
  VecSub(V[2],V[1],D);
  VecSub(V[3],V[1],E);
  VecSub(V[3],V[2],F);

  double Volume = VecDot( A, VecCross(B,C,Scratch) ) / 6.0;

  double TotalArea=0.0;
  TotalArea += VecNorm( VecCross(A,B,Scratch) );
  TotalArea += VecNorm( VecCross(A,C,Scratch) );
  TotalArea += VecNorm( VecCross(F,B,Scratch) );
  TotalArea += VecNorm( VecCross(F,D,Scratch) );
  TotalArea *= 0.5;

  double AvgLength = (   VecNorm(A) + VecNorm(B) 
                       + VecNorm(C) + VecNorm(D) 
                       + VecNorm(E) + VecNorm(F) 
                     ) / 6.0;

  if (VV) *VV=Volume;
  if (AA) *AA=TotalArea;
  if (LL) *LL=AvgLength;
  if (f)
   { 
     fprintf(f,"# Areas: (%.6f,%.6f,%.6f,%.6f)\n",
                0.5*VecNorm(VecCross(A,B,Scratch)),
                0.5*VecNorm(VecCross(A,C,Scratch)),
                0.5*VecNorm(VecCross(F,B,Scratch)),
                0.5*VecNorm(VecCross(F,D,Scratch)));
     fprintf(f,"# Lengths: (%.6f,%.6f,%.6f,%.6f,%6f,%6f)\n",

     fprintf(f,"%e %e %e \n",V[0][0],V[0][1],V[0][2]);
     fprintf(f,"%e %e %e \n",V[1][0],V[1][1],V[1][2]);
     fprintf(f,"%e %e %e \n",V[2][0],V[2][1],V[2][2]);
     fprintf(f,"%e %e %e \n",V[0][0],V[0][1],V[0][2]);
     fprintf(f,"\n");
     fprintf(f,"%e %e %e \n",V[0][0],V[0][1],V[0][2]);
     fprintf(f,"%e %e %e \n",V[1][0],V[1][1],V[1][2]);
     fprintf(f,"%e %e %e \n",V[3][0],V[3][1],V[3][2]);
     fprintf(f,"%e %e %e \n",V[0][0],V[0][1],V[0][2]);
     fprintf(f,"\n");
     fprintf(f,"%e %e %e \n",V[0][0],V[0][1],V[0][2]);
     fprintf(f,"%e %e %e \n",V[2][0],V[2][1],V[2][2]);
     fprintf(f,"%e %e %e \n",V[3][0],V[3][1],V[3][2]);
     fprintf(f,"%e %e %e \n",V[0][0],V[0][1],V[0][2]);
     fprintf(f,"\n");
     fprintf(f,"%e %e %e \n",V[1][0],V[1][1],V[1][2]);
     fprintf(f,"%e %e %e \n",V[2][0],V[2][1],V[2][2]);
     fprintf(f,"%e %e %e \n",V[3][0],V[3][1],V[3][2]);
     fprintf(f,"%e %e %e \n",V[1][0],V[1][1],V[1][2]);
     fprintf(f,"\n\n");

   };

  return Volume / (TotalArea*AvgLength*IDEAL_QUALITY_FACTOR);

}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *MeshFile=0;
  double QFThreshold=0.5;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"mesh",               PA_STRING, 1, 1, (void *)&MeshFile,    0, ".msh file"},
     {"meshfile",           PA_STRING, 1, 1, (void *)&MeshFile,    0, ".msh file"},
     {"QFThreshold",        PA_DOUBLE, 1, 1, (void *)&QFThreshold, 0, "quality-factor threshold"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (MeshFile==0)
   OSUsage(argv[0],OSArray,"--mesh / --meshfile option is mandatory");

  SWGVolume *V = new SWGVolume(MeshFile);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumVertices  = V->NumVertices;
  double *Vertices = (double *)memdup(V->Vertices, 3*NumVertices*sizeof(double));

  int NumTets     = 0;
  int *TetIndices = (int *)mallocEC( (16*V->NumTets)*sizeof(int) );

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double MaxQF=0.0, MinQF=1.0, AvgQF=0.0;
  int Subdivided=0;
  for(int nt=0; nt<V->NumTets; nt++)
   { 
     double QF = GetTetrahedronQualityFactor(V, nt);
     if (QF >= QFThreshold)
      { 
        TetIndices[ 4*NumTets + 0 ] = V->Tets[nt]->VI[0]; 
        TetIndices[ 4*NumTets + 1 ] = V->Tets[nt]->VI[1];
        TetIndices[ 4*NumTets + 2 ] = V->Tets[nt]->VI[2];
        TetIndices[ 4*NumTets + 3 ] = V->Tets[nt]->VI[3];
        NumTets++;
        MinQF  = fmin(MinQF, QF);
        MaxQF  = fmax(MaxQF, QF);
        AvgQF += QF;

      }
     else
      {
        Subdivided++;
        SWGTet *T=V->Tets[nt];
        printf("Tet %4i: QF=%.2e (subdividing...)\n",nt,QF);

        // add a new vertex corresponding to the centroid of this tet
        int nnv=NumVertices;
        NumVertices++;
        Vertices = (double *)realloc(Vertices, 3*NumVertices*sizeof(double));
        Vertices[3*nnv + 0] = T->Centroid[0];
        Vertices[3*nnv + 1] = T->Centroid[1];
        Vertices[3*nnv + 2] = T->Centroid[2];

        // add four new tets
        TetIndices[ 4*NumTets + 0 ] = nnv;
        TetIndices[ 4*NumTets + 1 ] = V->Tets[nt]->VI[1];
        TetIndices[ 4*NumTets + 2 ] = V->Tets[nt]->VI[2];
        TetIndices[ 4*NumTets + 3 ] = V->Tets[nt]->VI[3];
        NumTets++;

        TetIndices[ 4*NumTets + 0 ] = V->Tets[nt]->VI[0];
        TetIndices[ 4*NumTets + 1 ] = nnv;
        TetIndices[ 4*NumTets + 2 ] = V->Tets[nt]->VI[2];
        TetIndices[ 4*NumTets + 3 ] = V->Tets[nt]->VI[3];
        NumTets++;
 
        TetIndices[ 4*NumTets + 0 ] = V->Tets[nt]->VI[0]; 
        TetIndices[ 4*NumTets + 1 ] = V->Tets[nt]->VI[1];
        TetIndices[ 4*NumTets + 2 ] = nnv;
        TetIndices[ 4*NumTets + 3 ] = V->Tets[nt]->VI[3];
        NumTets++;
 
        TetIndices[ 4*NumTets + 0 ] = V->Tets[nt]->VI[0];
        TetIndices[ 4*NumTets + 1 ] = V->Tets[nt]->VI[1];
        TetIndices[ 4*NumTets + 2 ] = V->Tets[nt]->VI[2];
        TetIndices[ 4*NumTets + 3 ] = nnv;
        NumTets++;

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        FILE *f=vfopen("Tet%i.out","w",nt);
        double *Verts[4];

        Verts[0] = Vertices + 3*(V->Tets[nt]->VI[0]);
        Verts[1] = Vertices + 3*(V->Tets[nt]->VI[1]);
        Verts[2] = Vertices + 3*(V->Tets[nt]->VI[2]);
        Verts[3] = Vertices + 3*(V->Tets[nt]->VI[3]);
        double VV, AA, LL, QF2, VVTot=0.0;
        QF2 = GetTetQF(Verts, &VV, &AA, &LL, f);
        printf(" old: (V, A, L, QF)=(%.6f, %.6f, %.6f, %.6f)\n",VV,AA,LL,QF2);

        Verts[0] = Vertices + 3*nnv;
        Verts[1] = Vertices + 3*(V->Tets[nt]->VI[1]);
        Verts[2] = Vertices + 3*(V->Tets[nt]->VI[2]);
        Verts[3] = Vertices + 3*(V->Tets[nt]->VI[3]);
        QF2 = GetTetQF(Verts, &VV, &AA, &LL, f);
        VVTot += VV;
        printf(" new 1: (V, A, L, QF)=(%.6f, %.6f, %.6f, %.6f)\n",VV,AA,LL,QF2);

        Verts[0] = Vertices + 3*(V->Tets[nt]->VI[0]);
        Verts[1] = Vertices + 3*nnv;
        Verts[2] = Vertices + 3*(V->Tets[nt]->VI[2]);
        Verts[3] = Vertices + 3*(V->Tets[nt]->VI[3]);
        QF2 = GetTetQF(Verts, &VV, &AA, &LL, f);
        VVTot += VV;
        printf(" new 2: (V, A, L, QF)=(%.6f, %.6f, %.6f, %.6f)\n",VV,AA,LL,QF2);

        Verts[0] = Vertices + 3*(V->Tets[nt]->VI[0]);
        Verts[1] = Vertices + 3*(V->Tets[nt]->VI[1]);
        Verts[2] = Vertices + 3*nnv;
        Verts[3] = Vertices + 3*(V->Tets[nt]->VI[3]);
        QF2 = GetTetQF(Verts, &VV, &AA, &LL, f);
        VVTot += VV;
        printf(" new 3: (V, A, L, QF)=(%.6f, %.6f, %.6f, %.6f)\n",VV,AA,LL,QF2);

        Verts[0] = Vertices + 3*(V->Tets[nt]->VI[0]);
        Verts[1] = Vertices + 3*(V->Tets[nt]->VI[1]);
        Verts[2] = Vertices + 3*(V->Tets[nt]->VI[2]);
        Verts[3] = Vertices + 3*nnv;
        QF2 = GetTetQF(Verts, &VV, &AA, &LL, f);
        VVTot += VV;
        printf(" new 4: (V, A, L, QF)=(%.6f, %.6f, %.6f, %.6f)\n",VV,AA,LL,QF2);
        printf(" new total volume = %.6f\n",VVTot);
        printf(" \n");

      }; // if (QF >= QFThreshold ... else ... 

   }; // for(int nt=0; nt<V->NumTets; nt++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Subdivided==0)
   {
     printf("Subdivided 0 elements.\n");
   }
  else
   { 
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     FILE *f=vfopen("%s.ReMesh.msh","w",GetFileBase(MeshFile));
     fprintf(f,"$MeshFormat\n");
     fprintf(f,"2.2 0 8\n");
     fprintf(f,"$Nodes\n");
     fprintf(f,"%i\n",NumVertices);
     for(int nv=0; nv<NumVertices; nv++)
      fprintf(f,"%i %e %e %e \n",nv+1,Vertices[3*nv+0], Vertices[3*nv+1], Vertices[3*nv+2]);
     fprintf(f,"$EndNodes\n");
     fprintf(f,"$Elements\n");
     fprintf(f,"%i\n",NumTets);
     for(int nt=0; nt<NumTets; nt++)
      fprintf(f,"%i 4 2 0 1 %i %i %i %i\n",nt+1,
                 TetIndices[4*nt+0]+1, TetIndices[4*nt+1]+1,
                 TetIndices[4*nt+2]+1, TetIndices[4*nt+3]+1);
     fprintf(f,"$EndElements\n");
     fclose(f);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("Subdivided %i elements.\n",Subdivided);
     printf("New mesh file written to %s.remsh.\n",GetFileBase(MeshFile));
     
   };
  printf("Thank you for your support.\n");
 
}
