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
 * CreateBNEQData.cc -- a utility function to initialize a 
 *                   -- buff-neq data structure for a given
 *                   -- run of the code
 *
 * homer reid      -- 2/2012
 *
 */

#include "buff-neq.h"


/***************************************************************/
/***************************************************************/
/***************************************************************/
BNEQData *CreateBNEQData(char *GeoFile, char *TransFile, int QuantityFlags,
                         char *FileBase, 
                         bool DoOPFT, bool DoEMTPFT, bool DoMomentPFT,
                         int DSIPoints, double DSIRadius, char *DSIMesh, 
                         int DSIPoints2)
{

  BNEQData *BNEQD=(BNEQData *)mallocEC(sizeof(*BNEQD));

  /*--------------------------------------------------------------*/
  /*-- try to create the SWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  SWGGeometry *G=new SWGGeometry(GeoFile);
  BNEQD->G=G;
  
  if (FileBase)
   BNEQD->FileBase = strdup(FileBase);
  else
   FileBase = BNEQD->FileBase = strdup(GetFileBase(G->GeoFileName));

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  BNEQD->GTCList=ReadTransFile(TransFile, &(BNEQD->NumTransformations));
  char *ErrMsg=G->CheckGTCList(BNEQD->GTCList, BNEQD->NumTransformations);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  /*--------------------------------------------------------------*/
  /*- figure out which quantities were specified                 -*/
  /*--------------------------------------------------------------*/
  BNEQD->QuantityFlags=QuantityFlags;
  BNEQD->NQ=0;
  if ( QuantityFlags & QFLAG_PABS    ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_PRAD    ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_XFORCE  ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_YFORCE  ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_ZFORCE  ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_XTORQUE ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_YTORQUE ) BNEQD->NQ++;
  if ( QuantityFlags & QFLAG_ZTORQUE ) BNEQD->NQ++;
  BNEQD->NONQ = G->NumObjects * BNEQD->NQ; 
  BNEQD->NTNONQ = BNEQD->NumTransformations * BNEQD->NONQ;

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse -*/
  /*- chunks of the VIE matrices for multiple geometrical        -*/
  /*- transformations.                                           -*/
  /*--------------------------------------------------------------*/
  int NO=G->NumObjects;
  BNEQD->VBlocks = (SMatrix **)mallocEC(NO*sizeof(SMatrix *));
  BNEQD->RBlocks = (SMatrix **)mallocEC(NO*sizeof(SMatrix *));
  BNEQD->GBlocks = (HMatrix ***)mallocEC(NO*sizeof(HMatrix **));
  for(int no=0; no<NO; no++)
   { 
     int NBF  = G->Objects[no]->NumInteriorFaces;

     BNEQD->VBlocks[no]  = new SMatrix(NBF, NBF, LHM_COMPLEX);
     BNEQD->RBlocks[no]    = new SMatrix(NBF, NBF, LHM_REAL);

     BNEQD->VBlocks[no]->BeginAssembly(MAXOVERLAP);
     BNEQD->RBlocks[no]->BeginAssembly(MAXOVERLAP);

     BNEQD->GBlocks[no]  = (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
     int noMate=G->Mate[no];
     if (noMate!=-1)
      BNEQD->GBlocks[no][no] = BNEQD->GBlocks[noMate][noMate];
     else 
      BNEQD->GBlocks[no][no] = new HMatrix(NBF, NBF, LHM_COMPLEX);

     for(int nop=no+1; nop<NO; nop++)
      { int NBFP = G->Objects[nop]->NumInteriorFaces;
        BNEQD->GBlocks[no][nop] = new HMatrix(NBF, NBFP, LHM_COMPLEX);
      };
   };
   
  int NBFTot=G->TotalBFs;
  for(int n=0; n<3; n++)
   BNEQD->WorkMatrix[n] = new HMatrix(NBFTot, NBFTot, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*- compute and invert the overlap matrix (gram matrix)---------*/
  /*--------------------------------------------------------------*/
  HMatrix *SInverse=BNEQD->SInverse=new HMatrix(NBFTot, NBFTot, LHM_COMPLEX);
  SInverse->Zero();
  for(int no=0; no<NO; no++)
   { 
     SWGVolume *O = G->Objects[no];
     int Offset   = G->BFIndexOffset[no];
     for(int nfa=0; nfa<O->NumInteriorFaces; nfa++)
      { 
        int nfbList[MAXOVERLAP];
        double Entries[MAXOVERLAP];
        int NNZ=GetOverlapElements(O, nfa, nfbList, Entries);
        for(int nnz=0; nnz<NNZ; nnz++)
         SInverse->SetEntry(Offset+nfa, Offset+nfbList[nnz], Entries[nnz]);
      };
   };
  SInverse->LUFactorize();
  SInverse->LUInvert();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BNEQD->NumPFTMethods=0;
  int nPFT=0;
  if (DoOPFT)
   { BNEQD->PFTMethods[nPFT]   = SCUFF_PFT_OVERLAP;
     BNEQD->DSIPoints[nPFT]    = 0;
     BNEQD->PFTFileNames[nPFT] = vstrdup("%s.SIFlux.OPFT",FileBase);
     nPFT++;
   };
  if (DoMomentPFT)
   { BNEQD->PFTMethods[nPFT]   = SCUFF_PFT_MOMENTS;
     BNEQD->DSIPoints[nPFT]    = 0;
     BNEQD->PFTFileNames[nPFT] = vstrdup("%s.SIFlux.MomentPFT",FileBase);
     nPFT++;
   };
  if (DSIPoints>0 || DSIMesh)
   { BNEQD->PFTMethods[nPFT]   = SCUFF_PFT_DSI;
     BNEQD->DSIPoints[nPFT]    = DSIPoints;
     if (DSIPoints>0)
      BNEQD->PFTFileNames[nPFT] = vstrdup("%s.SIFlux.DSI%i",FileBase,DSIPoints);
     else
      BNEQD->PFTFileNames[nPFT] = vstrdup("%s.SIFlux.DSI.%s",FileBase,DSIMesh);
     nPFT++;
   };
  if (DSIPoints2>0)
   { BNEQD->PFTMethods[nPFT]   = SCUFF_PFT_DSI;
     BNEQD->DSIPoints[nPFT]    = DSIPoints2;
     BNEQD->PFTFileNames[nPFT] = vstrdup("%s.SIFlux.DSI%i",FileBase,DSIPoints2);
     nPFT++;
   };
  if (DoEMTPFT || nPFT==0)
   { BNEQD->PFTMethods[nPFT]   = SCUFF_PFT_EMT;
     BNEQD->DSIPoints[nPFT]    = 0;
     BNEQD->PFTFileNames[nPFT] = vstrdup("%s.SIFlux.EMTPFT",FileBase);
     nPFT++;
   };
  BNEQD->NumPFTMethods=nPFT;

  BNEQD->DSIOmegaPoints=0;

  BNEQD->pftOptions = (PFTOptions *)malloc(sizeof(PFTOptions));
  InitPFTOptions(BNEQD->pftOptions);
  BNEQD->pftOptions->DSIMesh=DSIMesh;
  BNEQD->pftOptions->DSIRadius=DSIRadius;
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BNEQD->TEnvironment=0.0;
  BNEQD->TemperatureSVTs=(SVTensor **)mallocEC(NO*sizeof(SVTensor *));
  BNEQD->TAvg=(double *)mallocEC(NO*sizeof(double));
  
  /*--------------------------------------------------------------*/
  /* write file preambles ----------------------------------------*/
  /*--------------------------------------------------------------*/
  time_t MyTime;
  struct tm *MyTm;
  char TimeString[30];
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  for(int n=0; n<BNEQD->NumPFTMethods; n++)
   { FILE *f=fopen(BNEQD->PFTFileNames[n],"a");
     fprintf(f,"# buff-neq run on %s (%s)\n",GetHostName(),TimeString);
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1 omega \n");
     fprintf(f,"# 2 transform tag\n");
     fprintf(f,"# 3 (sourceObject,destObject) \n");
     int nq=4;
     fprintf(f,"# %i absorbed power spectral density\n",nq++);
     fprintf(f,"# %i radiated power spectral density\n",nq++);
     fprintf(f,"# %i x-force flux spectral density\n",nq++);
     fprintf(f,"# %i y-force flux spectral density\n",nq++);
     fprintf(f,"# %i z-force flux spectral density\n",nq++);
     fprintf(f,"# %i x-torque flux spectral density\n",nq++);
     fprintf(f,"# %i y-torque flux spectral density\n",nq++);
     fprintf(f,"# %i z-torque flux spectral density\n",nq++);
     fclose(f);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("After CreateBNEQData: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  return BNEQD;

}
