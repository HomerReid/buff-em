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

void Invert3x3Matrix(cdouble M[3][3], cdouble W[3][3]);

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
        O->SVT->Evaluate(1.0, x0, Eps);

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
                                   const char *Tag,
                                   bool RealPart)
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
      // at the four tetrahedra vertices or at the 
      // tetrahedron centroid
      double EpsAvg[4];
      double *VV[4];
      bool UseCentroid=true;
      for(int n=0; n<4; n++)
       { 
         VV[n] = O->Vertices + 3*(T->VI[n]);

         double x0[3];
         if (UseCentroid)
          memcpy(x0, T->Centroid, 3*sizeof(double));
         else
          memcpy(x0, VV[n], 3*sizeof(double));
         if (O->GT)   O->GT->UnApply(x0);
         if (O->OTGT) O->OTGT->UnApply(x0);
 
         cdouble Eps[3][3];
         O->SVT->Evaluate(1.0, x0, Eps);

         if (RealPart)
          EpsAvg[n]=real(Eps[0][0] + Eps[1][1] + Eps[2][2])/3.0;
         else
          EpsAvg[n]=imag(Eps[0][0] + Eps[1][1] + Eps[2][2])/3.0;
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

void SWGGeometry::PlotPermittivity(const char *FileName,
                                   const char *Tag)
{ 
  (void)Tag;
  PlotPermittivity(FileName,"Re(Eps)",true);
  PlotPermittivity(FileName,"Im(Eps)",false);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool SWGGeometry::GetJE(int no, int nt, double *X0,
                        HVector *JVector, cdouble Omega,
                        cdouble JE[6])
{
  SWGVolume *O = Objects[no];
  int Offset   = BFIndexOffset[no];
  SWGTet *T    = O->Tets[nt];

  // sum contributions of up to 4 basis functions to get
  // current density at X0
  JE[0]=JE[1]=JE[2]=JE[3]=JE[4]=JE[5]=0.0;
  for(int iF=0; iF<4; iF++)
   { 
     int nf = T->FI[iF];
     if (nf >= O->NumInteriorFaces) continue;

     SWGFace *F = O->Faces[nf];
     cdouble JAlpha = JVector->GetEntry(Offset + nf);
     double PreFac = F->Area / (3.0*T->Volume);
     double Sign = (nt == F->iPTet) ? 1.0 : -1.0;
     int iQ      = (nt == F->iPTet) ? F->iQP : F->iQM;
     double *Q   = O->Vertices + 3*iQ;
     JE[0] += Sign*PreFac*JAlpha*(X0[0] - Q[0]);
     JE[1] += Sign*PreFac*JAlpha*(X0[1] - Q[1]);
     JE[2] += Sign*PreFac*JAlpha*(X0[2] - Q[2]);
   };

  // get total E-field from J
  // E=(iZ_0/k) * Chi^{-1} * J
  cdouble Chi[3][3], InvChi[3][3];
  O->SVT->Evaluate(Omega, X0, Chi);
  Chi[0][0] -= 1.0;
  Chi[1][1] -= 1.0;
  Chi[2][2] -= 1.0;
  Invert3x3Matrix(Chi, InvChi);

  cdouble PreFactor = II*ZVAC / Omega;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    JE[3+Mu] += PreFactor*InvChi[Mu][Nu]*JE[Nu];

  return true;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetTetVolume(double *V1, double *V2, double *V3, double *V4)
{
  // volume = (1/6) A \cdot (B x C)
  double A[3], B[3], C[3];
  for(int i=0; i<3; i++)
   { A[i] = V2[i] - V1[i];
     B[i] = V3[i] - V1[i];
     C[i] = V4[i] - V1[i];
   };
  return (  A[0] * (B[1]*C[2] - B[2]*C[1])
           +A[1] * (B[2]*C[0] - B[0]*C[2])
           +A[2] * (B[0]*C[1] - B[1]*C[0])
         ) / 6.0;
}

bool PointInTet(double *X0, SWGTet *T, double *Vertices)
{
  double *V[4];
  V[0] = Vertices + 3*T->VI[0];
  V[1] = Vertices + 3*T->VI[1];
  V[2] = Vertices + 3*T->VI[2];
  V[3] = Vertices + 3*T->VI[3];

  double VolumeSum=0.0;
  VolumeSum= fabs(GetTetVolume(X0,V[0],V[1],V[2]))
            +fabs(GetTetVolume(X0,V[0],V[2],V[3]))
            +fabs(GetTetVolume(X0,V[1],V[2],V[3]))
            +fabs(GetTetVolume(X0,V[0],V[1],V[3]));

  return ( (VolumeSum-T->Volume) < 1.0e-4 * T->Volume);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool SWGGeometry::GetJE(double *X0, HVector *JVector, cdouble Omega,
                        cdouble JE[6])
{
  for(int no=0; no<NumObjects; no++)
   for(int nf=0; nf<Objects[no]->NumInteriorFaces; nf++)
    { SWGVolume *O = Objects[no];
      SWGFace   *F = O->Faces[nf];
      double    R2 = F->Radius * F->Radius;
      if ( VecDistance2(X0, F->Centroid) > R2 )
       continue;

      if ( PointInTet(X0, O->Tets[F->iPTet], O->Vertices) )
       return GetJE(no, F->iPTet, X0, JVector, Omega, JE);
      if ( PointInTet(X0, O->Tets[F->iMTet], O->Vertices) )
       return GetJE(no, F->iMTet, X0, JVector, Omega, JE); 
    };
  memset(JE,0,6*sizeof(cdouble));
  return false;
}

/***************************************************************/
/* XJEMatrix[nt, 0,1,2] = coordinates of tet #n centroid       */
/* XJEMatrix[nt, 3,4,5] = {Jx,Jy,Jz} at tet #n centroid        */
/* XJEMatrix[nt, 6,7,8] = {Ex,Ey,Ez} at tet #n centroid        */
/***************************************************************/
HMatrix *SWGGeometry::GetXJEMatrix(HVector *JVector, cdouble Omega, HMatrix *XJEMatrix)
{
  int TotalTets=0;
  for(int no=0; no<NumObjects; no++)
   TotalTets+=Objects[no]->NumTets;

  /***************************************************************/
  /* (re)allocate matrix as necessary ****************************/
  /***************************************************************/
  if ( XJEMatrix &&
        (   XJEMatrix->NR!=TotalTets
         || XJEMatrix->NC!=9
         || XJEMatrix->RealComplex!=LHM_COMPLEX
        )
     )
   { Warn("wrong-size matrix in GetXJEMatrix (reallocating)");
     delete XJEMatrix;
     XJEMatrix=0;
   };
  if (XJEMatrix==0)
   XJEMatrix=new HMatrix(TotalTets, 9, LHM_COMPLEX);

  /***************************************************************/
  /* get current density and E-field at tetrahedron centroids    */
  /***************************************************************/
  for(int no=0, ntt=0; no<NumObjects; no++)
   for(int nt=0; nt<Objects[no]->NumTets; nt++, ntt++)
    { 
      double *X0 = Objects[no]->Tets[nt]->Centroid;
      cdouble JE[6];
      GetJE(no, nt, X0, JVector, Omega, JE);
      XJEMatrix->SetEntriesD(ntt, "0:2", X0);
      XJEMatrix->SetEntries(ntt,  "3:8", JE);
    };

  return XJEMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGGeometry::PlotCurrentDistribution(HVector *JVector, cdouble Omega,
                                          const char *format, ...)
{
  /***************************************************************/
  /***************************************************************/ 
  /***************************************************************/ 
  va_list ap;
  char FileBase[MAXSTR];
  va_start(ap,format);
  vsnprintfEC(FileBase,MAXSTR,format,ap);
  va_end(ap);
  char *p=strrchr(FileBase,'.');
  if (p && !strcasecmp(p,".pp"))
   *p=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  static HMatrix *XJEMatrix=0;
  XJEMatrix=GetXJEMatrix(JVector, Omega, XJEMatrix);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=vfopen("%s.XJE.dat","r",FileBase);
  if (!f)
   { f=vfopen("%s.XJE","w",FileBase);
     fprintf(f,"# 1,2,3,4  x,y,z,Omega\n");
     fprintf(f,"#  5,6     real,imag Jx\n");
     fprintf(f,"#  7,8     real,imag Jy\n");
     fprintf(f,"#  9,10    real,imag Jz\n");
     fprintf(f,"# 11,12    real,imag Ex\n");
     fprintf(f,"# 13,14    real,imag Ey\n");
     fprintf(f,"# 15,16    real,imag Ez\n");
     fprintf(f,"\n");
   };
  fclose(f);

  f=vfopen("%s.XJE.dat","a",FileBase);
  for(int nr=0; nr<XJEMatrix->NR; nr++)
   { double X0[3];
     cdouble JE[6]; 
     XJEMatrix->GetEntriesD(nr,"0:2",X0);
     XJEMatrix->GetEntries(nr,"3:8",JE);
     fprintVec(f,X0,3);
     fprintf(f,"%s ",z2s(Omega));
     fprintVecCR(f,JE,6);
   };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  f=vfopen("%s.XJE.pp","a",FileBase);
  if (!f) return;

  fprintf(f,"View \"Real J (Omega=%s)\" {\n",z2s(Omega));
  for(int nr=0; nr<XJEMatrix->NR; nr++)
   { 
     double X0[3];
     double J[3];
     XJEMatrix->GetEntriesD(nr,"0:2",X0);
     J[0] = real(XJEMatrix->GetEntry(nr,3));
     J[1] = real(XJEMatrix->GetEntry(nr,4));
     J[2] = real(XJEMatrix->GetEntry(nr,5));
     WriteVP(X0, J, f);
   };
  fprintf(f,"};\n");

  fprintf(f,"View \"Imag J (Omega=%s)\" {\n",z2s(Omega));
  for(int nr=0; nr<XJEMatrix->NR; nr++)
   { 
     double X0[3];
     double J[3];
     XJEMatrix->GetEntriesD(nr,"0:2",X0);
     J[0] = imag(XJEMatrix->GetEntry(nr,3));
     J[1] = imag(XJEMatrix->GetEntry(nr,4));
     J[2] = imag(XJEMatrix->GetEntry(nr,5));
     WriteVP(X0, J, f);
   };
  fprintf(f,"};\n");

  fprintf(f,"View \"Real E (Omega=%s)\" {\n",z2s(Omega));
  for(int nr=0; nr<XJEMatrix->NR; nr++)
   { 
     double X0[3];
     double E[3];
     XJEMatrix->GetEntriesD(nr,"0:2",X0);
     E[0] = real(XJEMatrix->GetEntry(nr,6));
     E[1] = real(XJEMatrix->GetEntry(nr,7));
     E[2] = real(XJEMatrix->GetEntry(nr,8));
     WriteVP(X0, E, f);
   };
  fprintf(f,"};\n");

  fprintf(f,"View \"Imag E (Omega=%s)\" {\n",z2s(Omega));
  for(int nr=0; nr<XJEMatrix->NR; nr++)
   { 
     double X0[3];
     double E[3];
     XJEMatrix->GetEntriesD(nr,"0:2",X0);
     E[0] = imag(XJEMatrix->GetEntry(nr,6));
     E[1] = imag(XJEMatrix->GetEntry(nr,7));
     E[2] = imag(XJEMatrix->GetEntry(nr,8));
     WriteVP(X0, E, f);
   };
  fprintf(f,"};\n");

  fclose(f);
}

} // namespace buff
