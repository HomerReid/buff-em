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
 * OutputModules.cc -- various types of 'output modules' for buff-scatter 
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "buff-scatter.h"

#define ABSTOL   0.0
#define RELTOL   5.0e-2
#define MAXEVALS 20000

#define II cdouble(0.0,1.0)

using namespace buff;
using namespace scuff;

namespace buff {

void GetMoments(SWGGeometry *G, int no, cdouble Omega,
                HVector *JVector,
                cdouble p[3], cdouble m[3], cdouble Qp[3]);

void AddIFContributionsToEMTPFT(SWGGeometry *G, HVector *JVector,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix);

HMatrix *GetEMTPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HVector *RHSVector,
                   HMatrix *DMatrix, HMatrix *PFTMatrix);

void GetMomentPFT(SWGGeometry *G, int no, cdouble Omega,
                  IncField *IF, HVector *JVector, HMatrix *DMatrix,
                  HMatrix *PFTMatrix, bool KeepQpTerm, 
                  char *FileBase);

void GetMomentPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                  HVector *JVector, HMatrix *DMatrix,
                  HMatrix *PFTMatrix, bool KeepQpTerm, 
                  char *FileBase);
               }

/***************************************************************/
/* compute scattered and total fields at a user-specified list */
/* of evaluation points                                        */
/***************************************************************/
void ProcessEPFile(BSData *BSD, char *EPFileName)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SWGGeometry *G  = BSD->G;
  IncField *IF    = BSD->IF;
  HVector  *J     = BSD->J; 
  cdouble  Omega  = BSD->Omega;

  /*--------------------------------------------------------------*/
  /*- try to read eval points from file --------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(EPFileName,LHM_TEXT,"-ncol 3");
  if (XMatrix->ErrMsg)
   { fprintf(stderr,"Error processing EP file: %s\n",XMatrix->ErrMsg);
     delete XMatrix;
     return;
   };

  /*--------------------------------------------------------------*/
  /*- get components of scattered and incident fields            -*/
  /*--------------------------------------------------------------*/
  Log("Evaluating fields at points in file %s...",EPFileName);

  HMatrix *FMatrix1 = G->GetFields( 0, J, Omega, XMatrix); // scattered
  HMatrix *FMatrix2 = G->GetFields(IF, 0, Omega, XMatrix); // incident

  /*--------------------------------------------------------------*/
  /*- create .scattered and .total output files and write fields -*/
  /*--------------------------------------------------------------*/
  char buffer[MAXSTR];
  snprintf(buffer,MAXSTR,"%s.scattered",GetFileBase(EPFileName));
  FILE *f1=CreateUniqueFile(buffer,1);
  snprintf(buffer,MAXSTR,"%s.total",GetFileBase(EPFileName));
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
  delete XMatrix;
  delete FMatrix1;
  delete FMatrix2;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePFTFile(BSData *BSD, char *PFTFile, PFTOptions *Options,
                  int PFTMethod)
{ 
  /***************************************************************/
  /* write file preamble as necessary ****************************/
  /***************************************************************/
  FILE *f=fopen(PFTFile,"r");
  if (!f)
   { 
     f=fopen(PFTFile,"w");
     fprintf(f,"# data columns:    \n");
     fprintf(f,"# 1 frequency      \n");
     fprintf(f,"# 2 object label   \n");
     fprintf(f,"# 3 absorbed power (watts) \n");
     fprintf(f,"# 4 radiated power (watts) \n");
     fprintf(f,"# 5,6,7  x, y, z force (nanoNewtons)\n");
     fprintf(f,"# 8,9,10 x, y, z torque (nanoNewtons microns)\n");
     fprintf(f,"\n");
   };
  fclose(f);

  /***************************************************************/
  /* do the PFT calculation **************************************/
  /***************************************************************/
  SWGGeometry *G     = BSD->G;
  HVector *J         = BSD->J;
  cdouble Omega      = BSD->Omega;
  Options->IF        = BSD->IF;
  Options->RHSVector = BSD->RHS;
  Options->PFTMethod = PFTMethod;

  HMatrix *PFTMatrix =G->GetPFT(J, Omega, 0, Options);

  /***************************************************************/
  /* write results to output file ********************************/
  /***************************************************************/
  f=fopen(PFTFile,"a");
  for(int no=0; no<G->NumObjects; no++)
   { fprintf(f,"%e %s ",real(BSD->Omega),G->Objects[no]->Label);
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFTMatrix->GetEntryD(no,nq));
     fprintf(f,"\n");
   };
  fclose(f);

  delete PFTMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteEMTMomentPFTFile(BSData *BSD, char *PFTFile)
{ 
  /***************************************************************/
  /* write file preamble as necessary ****************************/
  /***************************************************************/
  FILE *f=fopen(PFTFile,"r");
  if (!f)
   { 
     f=fopen(PFTFile,"w");
     fprintf(f,"# data columns:          \n");
     fprintf(f,"# 1 frequency            \n");
     fprintf(f,"# 2 object label         \n");
     fprintf(f,"# 3-10  J--I PFT (EMT)   \n");
     fprintf(f,"# 11-18 J--J PFT (EMT)   \n");
     fprintf(f,"# 19-26 J--I PFT (moments) \n");
     fprintf(f,"# 27-34 J--J PFT (moments wo) \n");
     fprintf(f,"# 35-42 J--J PFT (moments w) \n");
     fprintf(f,"# 43-45 re (Px, Py, Pz)  \n");
     fprintf(f,"# 46-48 im (Px, Py, Pz)  \n");
     fprintf(f,"# 49-51 re (Mx, My, Mz)  \n");
     fprintf(f,"# 52-54 im (Mx, My, Mz)  \n");
     fprintf(f,"# 55 J--I time (EMT)\n");
     fprintf(f,"# 56 J--J time (EMT)\n");
     fprintf(f,"# 57 J--I time (moments)\n");
     fprintf(f,"# 58 J--J time (moments)\n");
     fprintf(f,"\n");
   };
  fclose(f);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G     = BSD->G;
  HVector *J         = BSD->J;
  cdouble Omega      = BSD->Omega;
  IncField *IF       = BSD->IF;
  HVector *RHSVector = BSD->RHS;

  /*--------------------------------------------------------------*/
  /*- get J-I and J-J contributions to EMTPFT --------------------*/
  /*--------------------------------------------------------------*/
  int NO=G->NumObjects;
  HMatrix *PFTMatrix1 = new HMatrix(NO, NUMPFT);
  HMatrix *PFTMatrix2 = new HMatrix(NO, NUMPFT);
  HMatrix *PFTMatrix3 = new HMatrix(NO, NUMPFT);
  HMatrix *PFTMatrix4 = new HMatrix(NO, NUMPFT);
  HMatrix *PFTMatrix5 = new HMatrix(NO, NUMPFT);
  PFTMatrix1->Zero();
  PFTMatrix2->Zero();
  PFTMatrix3->Zero();
  PFTMatrix4->Zero();
  PFTMatrix5->Zero();

  Tic();
  AddIFContributionsToEMTPFT(G, J, IF, Omega, PFTMatrix1);
  GetEMTPFT(G, Omega, 0, J, RHSVector, 0, PFTMatrix2);
  double EMTTime=Toc();

  Tic();
  GetMomentPFT(G, Omega, IF, J, 0, PFTMatrix3, false, 0);
  double MomentTime=Toc();

  GetMomentPFT(G, Omega,  0, J, 0, PFTMatrix4, false, 0);
  GetMomentPFT(G, Omega,  0, J, 0, PFTMatrix5,  true, 0);

  for(int no=0; no<NO; no++)
   for(int nq=0; nq<NUMPFT; nq++)
    PFTMatrix3->AddEntry(no,nq,-1.0*PFTMatrix4->GetEntry(no,nq));

  /***************************************************************/
  /* write results to output file ********************************/
  /***************************************************************/
  f=fopen(PFTFile,"a");
  for(int no=0; no<G->NumObjects; no++)
   { 
     fprintf(f,"%e %s ",real(BSD->Omega),G->Objects[no]->Label);

     char *FileBase = GetFileBase(G->GeoFileName);
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFTMatrix1->GetEntryD(no,nq));
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFTMatrix2->GetEntryD(no,nq));
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFTMatrix3->GetEntryD(no,nq));
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFTMatrix4->GetEntryD(no,nq));
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFTMatrix5->GetEntryD(no,nq));
      
     cdouble p[3], m[3], Qp[3];
     GetMoments(G, no, Omega, J, p, m, Qp);
     for(int Mu=0; Mu<3; Mu++) fprintf(f,"%e ",real(p[Mu]));
     for(int Mu=0; Mu<3; Mu++) fprintf(f,"%e ",imag(p[Mu]));
     for(int Mu=0; Mu<3; Mu++) fprintf(f,"%e ",real(m[Mu]));
     for(int Mu=0; Mu<3; Mu++) fprintf(f,"%e ",imag(m[Mu]));

     fprintf(f,"%e ",EMTTime);
     fprintf(f,"%e ",MomentTime);

     fprintf(f,"\n");
   };
  fclose(f);

  delete PFTMatrix1;
  delete PFTMatrix2;
  delete PFTMatrix3;
  delete PFTMatrix4;
  delete PFTMatrix5;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void WriteMomentFile(BSData *BSD, char *FileName)
{
  SWGGeometry *G    = BSD->G;
  HVector *JVector  = BSD->J;
  cdouble Omega     = BSD->Omega;

  FILE *f=fopen(FileName,"a");
  for(int no=0; no<G->NumObjects; no++)
   { 
     SWGVolume *O = G->Objects[no];
     fprintf(f,"%s %s ",z2s(Omega),O->Label);

     cdouble M[3]={0.0, 0.0, 0.0};
     int NBF    = G->Objects[no]->NumInteriorFaces;
     int Offset = G->BFIndexOffset[no];
     for(int nbf=0; nbf<NBF; nbf++)
      { 
        SWGFace *F = O->Faces[nbf];
        double *QP = O->Vertices + 3*(F->iQP);
        double *QM = O->Vertices + 3*(F->iQM);
        double A   = F->Area;

        cdouble PreFactor = A*(JVector->GetEntry(Offset + nbf)) / (II*Omega*4.0);

        M[0] += PreFactor*(QP[0] - QM[0]);
        M[1] += PreFactor*(QP[1] - QM[1]);
        M[2] += PreFactor*(QP[2] - QM[2]);
      };

     fprintf(f,"%e %e %e %e %e %e \n",
                real(M[0]),imag(M[0]),
                real(M[1]),imag(M[1]),
                real(M[2]),imag(M[2]));
   };
  fclose(f);
  
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteMomentFile(BSData *BSD, char *FileName)
{
  SWGGeometry *G    = BSD->G;
  HVector *JVector  = BSD->J;
  cdouble Omega     = BSD->Omega;

  FILE *f=fopen(FileName,"a");
  for(int no=0; no<G->NumObjects; no++)
   { 
     cdouble p[3], m[3], Qp[3];
     GetMoments(G, no, Omega, JVector, p, m, Qp);

     fprintf(f,"%s %s ",z2s(Omega),G->Objects[no]->Label);
     fprintf(f,"%e %e ",real(p[0]),imag(p[0]));
     fprintf(f,"%e %e ",real(p[1]),imag(p[1]));
     fprintf(f,"%e %e ",real(p[2]),imag(p[2]));
     fprintf(f,"%e %e ",real(m[0]),imag(m[0]));
     fprintf(f,"%e %e ",real(m[1]),imag(m[1]));
     fprintf(f,"%e %e ",real(m[2]),imag(m[2]));
     fprintf(f,"\n");
   };
  fclose(f);
  
}
