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
  char *FileBase  = BSD->FileBase;

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

  HMatrix *SFMatrix = G->GetFields( 0, J, Omega, XMatrix); // scattered
  HMatrix *IFMatrix = G->GetFields(IF, 0, Omega, XMatrix); // incident

  /*--------------------------------------------------------------*/
  /*- create .scattered and .total output files and write fields -*/
  /*--------------------------------------------------------------*/
  SetDefaultCD2SFormat("%+.8e %+.8e ");
  char OmegaStr[100];
  snprintf(OmegaStr,100,"%s",z2s(Omega));
  char *TransformLabel=BSD->TransformLabel;
  char *IFLabel=BSD->IFLabel;
  const char *Ext[2]={"scattered","total"};
  for(int ST=0; ST<2; ST++)
   { char OutFileName[MAXSTR];
     snprintf(OutFileName,MAXSTR,"%s.%s.%s",FileBase,GetFileBase(EPFileName),Ext[ST]);
     FILE *f=fopen(OutFileName,"a");
     fprintf(f,"# buff-scatter run on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# columns: \n");
     fprintf(f,"# 1,2,3   x,y,z (evaluation point coordinates)\n");
     fprintf(f,"# 4       omega (angular frequency)\n");
     int nc=5;
     if (TransformLabel)
      fprintf(f,"# %i       geometrical transform\n",nc++);
     if (IFLabel)
      fprintf(f,"# %i       incident field\n",nc++);
     fprintf(f,"# %02i,%02i   real, imag Ex\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Ex\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Ey\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Ez\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Hx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Hy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Hz\n",nc,nc+1); nc+=2;
     for(int nr=0; nr<SFMatrix->NR; nr++)
      { double X[3];
        cdouble EH[6];
        XMatrix->GetEntriesD(nr,":",X);
        SFMatrix->GetEntries(nr,":",EH);
        if (ST==1) 
         for(nc=0; nc<6; nc++) 
          EH[nc]+=IFMatrix->GetEntry(nr,nc);
        fprintf(f,"%+.8e %+.8e %+.8e ",X[0],X[1],X[2]);
        fprintf(f,"%s ",OmegaStr);
        if (TransformLabel) fprintf(f,"%s ",TransformLabel);
        if (IFLabel) fprintf(f,"%s ",IFLabel);
        fprintf(f,"%s %s %s   ",CD2S(EH[0]),CD2S(EH[1]),CD2S(EH[2]));
        fprintf(f,"%s %s %s\n", CD2S(EH[3]),CD2S(EH[4]),CD2S(EH[5]));
      };
     fclose(f);
   };

  delete XMatrix;
  delete SFMatrix;
  delete IFMatrix;

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
  char *TransformLabel = BSD->TransformLabel;
  char *IFLabel        = BSD->IFLabel;
  FILE *f=fopen(PFTFile,"r");
  if (!f)
   {  f=fopen(PFTFile,"w");
      fprintf(f,"# buff-scatter run on %s (%s)\n",GetHostName(),GetTimeString());
      fprintf(f,"# data file columns: \n");
      fprintf(f,"# 1   omega           (rad/sec) \n");
      int nc=2;
      if (TransformLabel)
       fprintf(f,"# %i   geometrical transform\n",nc++);
      if (IFLabel)
       fprintf(f,"# %i   incident field\n",nc++);
      fprintf(f,"#%2i   surface label \n",nc++);             
      fprintf(f,"#%2i   absorbed power  (watts)\n",nc++);
      fprintf(f,"#%2i   scattered power (watts)\n",nc++);
      fprintf(f,"#%2i   x-force         (nanonewtons)\n",nc++);
      fprintf(f,"#%2i   y-force         (nanonewtons)\n",nc++);
      fprintf(f,"#%2i   z-force         (nanonewtons)\n",nc++);
      fprintf(f,"#%2i   x-torque        (nanonewtons * microns)\n",nc++);
      fprintf(f,"#%2i   y-torque        (nanonewtons * microns)\n",nc++);
      fprintf(f,"#%2i   z-torque        (nanonewtons * microns)\n",nc++);
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

  HMatrix *PFTMatrix =G->GetPFTMatrix(J, Omega, Options);

  /***************************************************************/
  /* write results to output file ********************************/
  /***************************************************************/
  f=fopen(PFTFile,"a");
  for(int no=0; no<G->NumObjects; no++)
   { fprintf(f,"%e ",real(BSD->Omega));
     if (TransformLabel) 
      fprintf(f,"%s ",TransformLabel);
     if (IFLabel) 
      fprintf(f,"%s ",IFLabel);
     fprintf(f,"%s ",G->Objects[no]->Label);
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
void WriteMomentFile(BSData *BSD, char *FileName)
{
  char *TransformLabel = BSD->TransformLabel;
  char *IFLabel        = BSD->IFLabel;

  /***************************************************************/
  /* write file preamble on initial file creation ****************/
  /***************************************************************/
  FILE *f=fopen(FileName,"r");
  if (!f)
   { f=fopen(FileName,"w");
     if ( !f ) ErrExit("could not open file %s",FileName);
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1     angular frequency (3e14 rad/sec)\n");
     int nc=2;
     if (TransformLabel)
      fprintf(f,"# %i     geometrical transform\n",nc++);
     if (IFLabel)
      fprintf(f,"# %i     incident field\n",nc++);
     fprintf(f,"# %i     surface label\n",nc++);
     fprintf(f,"# %02i,%02i real,imag px (electric dipole moment)\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag py \n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag pz \n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag mx (magnetic dipole moment)\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag my \n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag mz \n",nc,nc+1); nc+=2;
     fprintf(f,"#\n");
   };
  fclose(f);

  SWGGeometry *G       = BSD->G;
  HVector *JVector     = BSD->J;
  cdouble Omega        = BSD->Omega;

  f=fopen(FileName,"a");
  for(int no=0; no<G->NumObjects; no++)
   { 
     cdouble p[3], m[3], Qp[3];
     GetMoments(G, no, Omega, JVector, p, m, Qp);

     fprintf(f,"%s ",z2s(Omega));
     if (TransformLabel)
      fprintf(f,"%s ",TransformLabel);
     if (IFLabel)
      fprintf(f,"%s ",IFLabel);
     fprintf(f,"%s ",G->Objects[no]->Label);
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
