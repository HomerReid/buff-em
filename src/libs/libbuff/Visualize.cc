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
 * Visualize.cc
 * 
 * homer reid -- 11/2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;
namespace buff {

#define II cdouble(0,1)

#define MAXSTR 1000

/************************************************************/
/* subroutines for emitting GMSH postprocessing code        */
/************************************************************/
/* vector point (otherwise known as 'arrow') */
void WriteVP(double *X, double *V, FILE *f)
{
  fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",X[0],X[1],X[2],V[0],V[1],V[2]);
}

/* scalar triangle */
void WriteST(double **VV, double Val, FILE *f)
{
  fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
             VV[0][0], VV[0][1], VV[0][2],
             VV[1][0], VV[1][1], VV[1][2],
             VV[2][0], VV[2][1], VV[2][2], 
             Val,Val,Val);
}

/* scalar triangle */
void WriteST(double **VV, cdouble Val, FILE *f)
{
  fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
             VV[0][0], VV[0][1], VV[0][2],
             VV[1][0], VV[1][1], VV[1][2],
             VV[2][0], VV[2][1], VV[2][2], 
             real(Val),real(Val),real(Val));
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGVolume::DrawBF(int nf, FILE *f)
{
  /***************************************************************/
  /***************************************************************/ 
  /***************************************************************/ 
  fprintf(f,"View \"BF%i\" {\n",nf);

  SWGTet *TP = Tets[ Faces[nf]->iPTet ];
  fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
             TP->Centroid[0],TP->Centroid[1],TP->Centroid[2],TP->Index);

  for(int nfP=0; nfP<4; nfP++)
   { 
     SWGFace *F = Faces[ TP->FI[nfP] ];
     double *V1 = Vertices + 3*(F->iV1);
     double *V2 = Vertices + 3*(F->iV2);
     double *V3 = Vertices + 3*(F->iV3);
   
     fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                V1[0],V1[1],V1[2], V2[0],V2[1],V2[2], V3[0],V3[1],V3[2],
                +1.0, +1.0, +1.0);
     fprintf(f,"T3(%e,%e,%e,0.0) {\"%i(%i)\"};\n",
                F->Centroid[0], F->Centroid[1], F->Centroid[2], TP->Index, nfP);
   };

  SWGTet *TM = Tets[ Faces[nf]->iMTet ];
  fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
             TM->Centroid[0],TM->Centroid[1],TM->Centroid[2],TM->Index);
  for(int nfM=0; nfM<4; nfM++)
   { 
     if ( TM->FI[nfM] == nf )
      continue;

     SWGFace *F = Faces[ TM->FI[nfM] ];
     double *V1 = Vertices + 3*(F->iV1);
     double *V2 = Vertices + 3*(F->iV2);
     double *V3 = Vertices + 3*(F->iV3);
   
     fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                V1[0],V1[1],V1[2], V2[0],V2[1],V2[2], V3[0],V3[1],V3[2],
                -1.0, -1.0, -1.0);
     fprintf(f,"T3(%e,%e,%e,0.0) {\"%i(%i)\"};\n",
                F->Centroid[0], F->Centroid[1], F->Centroid[2], 
                TM->Index, nfM);
   };
  fprintf(f,"};\n");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGVolume::DrawTet(int nt, FILE *f)
{
  /***************************************************************/
  /***************************************************************/ 
  /***************************************************************/ 
  fprintf(f,"View \"Tet%i\" {\n",nt);

  SWGTet *T = Tets[ nt ];
  fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
             T->Centroid[0],T->Centroid[1],T->Centroid[2],T->Index);

  for(int nf=0; nf<4; nf++)
   { 
     SWGFace *F = Faces[ T->FI[nf] ];
     double *V1 = Vertices + 3*(F->iV1);
     double *V2 = Vertices + 3*(F->iV2);
     double *V3 = Vertices + 3*(F->iV3);
   
     fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                V1[0],V1[1],V1[2], V2[0],V2[1],V2[2], V3[0],V3[1],V3[2],
                +1.0, +1.0, +1.0);
     fprintf(f,"T3(%e,%e,%e,0.0) {\"%i(%i)\"};\n",
                F->Centroid[0], F->Centroid[1], F->Centroid[2], T->Index, nf);
   };
  fprintf(f,"};\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGVolume::DrawFace(int nf, FILE *f)
{

  SWGFace *F = Faces[nf];
  if (F->Index != nf)
   printf("Bawonkatage foryaf! %i \\ne %i \n",nf,F->Index);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"Face%i\" {\n",nf);

  fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
             F->Centroid[0],F->Centroid[1],F->Centroid[2],F->Index);

  double *V1 = Vertices + 3*(F->iV1);
  double *V2 = Vertices + 3*(F->iV2);
  double *V3 = Vertices + 3*(F->iV3);
   
  fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
             V1[0],V1[1],V1[2], V2[0],V2[1],V2[2], V3[0],V3[1],V3[2],
             +1.0, +1.0, +1.0);
  fprintf(f,"};\n");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::WritePPMesh(const char *FileName, 
                              const char *Tag)
{
  /***************************************************************/
  /***************************************************************/ 
  /***************************************************************/ 
  FILE *f=fopen(FileName,"a");
  if (!f) return;
  fprintf(f,"View \"%s\" {\n",Tag);

  /***************************************************************/
  /* draw all exterior faces on all objects, with each face      */
  /* colored by the average of the diagonal elements of the      */
  /* permittivity tensor at its centroid, evaluated at Omega=1.0 */
  /***************************************************************/
  SWGVolume *O; 
  int nf;
  for(int no=0; no<NumObjects; no++)
   for(O=Objects[no], nf=O->NumInteriorFaces; nf<O->NumTotalFaces; nf++)
   { 
     SWGFace *F = O->Faces[nf];
     double *VV[3];
     VV[0] = O->Vertices + 3*(F->iV1);
     VV[1] = O->Vertices + 3*(F->iV2);
     VV[2] = O->Vertices + 3*(F->iV3);

     // get average of diagonal permittivity elements
     // at the three face vertices
     double EpsAvg[3];
     for(int n=0; n<3; n++)
      { 
        double x0[3];
        memcpy(x0, VV[n], 3*sizeof(double));
        if (O->GT)   O->GT->UnApply(x0);
        if (O->OTGT) O->OTGT->UnApply(x0);

        cdouble Eps[3][3];
        O->MP->GetEps(1.0, x0, Eps);

        EpsAvg[n]=real(Eps[0][0] + Eps[1][1] + Eps[2][2])/3.0;
      };
   
     fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                VV[0][0],VV[0][1],VV[0][2],
                VV[1][0],VV[1][1],VV[1][2],
                VV[2][0],VV[2][1],VV[2][2],
                EpsAvg[0], EpsAvg[1], EpsAvg[2]);
   };
  fprintf(f,"};\n");
  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::PlotPermittivity(const char *FileName,
                                   const char *Tag)
{
  FILE *f=fopen(FileName,"a");
  if (!f) return;
  fprintf(f,"View \"%s\" {\n",Tag);
  for(int no=0; no<NumObjects; no++)
   for(int nt=0; nt<Objects[no]->NumTets; nt++)
    { 
      SWGVolume *O = Objects[no];
      SWGTet *T = O->Tets[nt];

      // get average of diagonal permittivity elements
      // at the four tetrahedra vertices
      double EpsAvg[4];
      double *VV[4];
      for(int n=0; n<4; n++)
       { 
         VV[n] = O->Vertices + 3*(T->VI[n]);

         double x0[3];
         memcpy(x0, VV[n], 3*sizeof(double));
         if (O->GT)   O->GT->UnApply(x0);
         if (O->OTGT) O->OTGT->UnApply(x0);
 
         cdouble Eps[3][3];
         O->MP->GetEps(1.0, x0, Eps);

         EpsAvg[n]=real(Eps[0][0] + Eps[1][1] + Eps[2][2])/3.0;
       };
   
     fprintf(f,"SS(%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e,%e};\n",
                VV[0][0],VV[0][1],VV[0][2],
                VV[1][0],VV[1][1],VV[1][2],
                VV[2][0],VV[2][1],VV[2][2],
                VV[3][0],VV[3][1],VV[3][2],
                EpsAvg[0], EpsAvg[1], EpsAvg[2], EpsAvg[3]);
    };
  fprintf(f,"};\n");
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::PlotCurrentDistribution(const char *FileName,
                                          HVector *J, 
                                          const char *Tag, ...)
{
  /***************************************************************/
  /***************************************************************/ 
  /***************************************************************/ 
  FILE *f=fopen(FileName,"a");
  if (!f) return;

  va_list ap;
  char buffer[MAXSTR];
  va_start(ap,Tag);
  vsnprintfEC(buffer,MAXSTR,Tag,ap);
  va_end(ap);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"Real(%s)\" {\n",buffer);
  for(int no=0; no<NumObjects; no++)
   for(int nt=0; nt<Objects[no]->NumTets; nt++)
    { 
      SWGVolume *O = Objects[no];
      int Offset   = BFIndexOffset[no];
      SWGTet *T    = O->Tets[nt];

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      double JCentroid[3]={0.0, 0.0, 0.0};
      for(int iF=0; iF<4; iF++)
       { 
         int nf = T->FI[iF];
         if (nf >= O->NumInteriorFaces) continue;

         SWGFace *F = O->Faces[nf];
         double JAlpha = real(J->GetEntry(Offset + nf));
         double PreFac = F->Area / (3.0*T->Volume);
         double Sign=1.0; 
         double *Q=0;
         if ( F->iPTet == nt)
          { Sign = +1.0;
            Q    = O->Vertices + 3*F->iQP;
          }
         else if ( F->iMTet == nt)
          { Sign = -1.0;
            Q    = O->Vertices + 3*F->iQM;
          }
         else 
          ErrExit("%s:%i: internal error",__FILE__,__LINE__);
   
         JCentroid[0] += Sign*PreFac*JAlpha*(T->Centroid[0] - Q[0]);
         JCentroid[1] += Sign*PreFac*JAlpha*(T->Centroid[1] - Q[1]);
         JCentroid[2] += Sign*PreFac*JAlpha*(T->Centroid[2] - Q[2]);
       };
      WriteVP(T->Centroid, JCentroid, f);
    };
  fprintf(f,"};\n");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"Imag(%s)\" {\n",buffer);
  for(int no=0; no<NumObjects; no++)
   for(int nt=0; nt<Objects[no]->NumTets; nt++)
    { 
      SWGVolume *O = Objects[no];
      int Offset   = BFIndexOffset[no];
      SWGTet *T    = O->Tets[nt];

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      double JCentroid[3]={0.0, 0.0, 0.0};
      for(int iF=0; iF<4; iF++)
       { 
         int nf = T->FI[iF];
         if (nf >= O->NumInteriorFaces) continue;

         SWGFace *F = O->Faces[nf];
         double JAlpha = imag(J->GetEntry(Offset + nf));
         double PreFac = F->Area / (3.0*T->Volume);
         double Sign=1.0;
         double *Q=0;
         if ( F->iPTet == nt)
          { Sign = +1.0;
            Q    = O->Vertices + 3*F->iQP;
          }
         else if ( F->iMTet == nt)
          { Sign = -1.0;
            Q    = O->Vertices + 3*F->iQM;
          }
         else 
          ErrExit("%s:%i: internal error",__FILE__,__LINE__);
   
         JCentroid[0] += Sign*PreFac*JAlpha*(T->Centroid[0] - Q[0]);
         JCentroid[1] += Sign*PreFac*JAlpha*(T->Centroid[1] - Q[1]);
         JCentroid[2] += Sign*PreFac*JAlpha*(T->Centroid[2] - Q[2]);
       };
      WriteVP(T->Centroid, JCentroid, f);
    };
  fprintf(f,"};\n");

  fclose(f);
}

} // namespace buff
