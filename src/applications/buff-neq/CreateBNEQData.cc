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
BNEQData *CreateBNEQData(char *GeoFile, char *TransFile,
                         char *TemperatureFile, int QuantityFlags,
                         char *pFileBase)
{

  BNEQData *BNEQD=(BNEQData *)mallocEC(sizeof(*BNEQD));

  /*--------------------------------------------------------------*/
  /*-- try to create the SWGGeometry -----------------------------*/
  /*--------------------------------------------------------------*/
  SWGGeometry *G=new SWGGeometry(GeoFile);
  BNEQD->G=G;
  
  if (pFileBase)
   BNEQD->FileBase = strdup(pFileBase);
  else
   BNEQD->FileBase = strdup(GetFileBase(G->GeoFileName));

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
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (TemperatureFile)
   BNEQD->Temperature = new IHAIMatProp(TemperatureFile);
  else
   BNEQD->Temperature = 0;

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
  BNEQD->VInv    = (SMatrix **)mallocEC(NO*sizeof(SMatrix *));
  BNEQD->Sigma   = (SMatrix **)mallocEC(NO*sizeof(SMatrix *));
  BNEQD->GBlocks = (HMatrix ***)mallocEC(NO*sizeof(HMatrix **));
  for(int no=0; no<NO; no++)
   { 
     int NBF  = G->Objects[no]->NumInteriorFaces;
     BNEQD->VInv[no]  = new SMatrix(NBF, NBF, LHM_COMPLEX);
     BNEQD->Sigma[no] = new SMatrix(NBF, NBF, LHM_REAL);
     BNEQD->GBlocks[no] = (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
     for(int nop=no; nop<NO; nop++)
      { int NBFP = G->Objects[nop]->NumInteriorFaces;
        BNEQD->GBlocks[no][nop] = new HMatrix(NBF, NBFP, LHM_COMPLEX);
      };
   };
  BNEQD->SMatricesInitialized=false;
   
  int NBFTot=G->TotalBFs;
  for(int n=0; n<3; n++)
   BNEQD->WorkMatrix[n] = new HMatrix(NBFTot, NBFTot, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /* write file preambles ----------------------------------------*/
  /*--------------------------------------------------------------*/
  time_t MyTime;
  struct tm *MyTm;
  char TimeString[30];
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  FILE *f=vfopen("%s.SRFlux","a",BNEQD->FileBase);
  fprintf(f,"\n");
  fprintf(f,"# buff-neq run on %s (%s)\n",GetHostName(),TimeString);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 omega \n");
  fprintf(f,"# 2 transform tag\n");
  fprintf(f,"# 3 (sourceObject,destObject) \n");
  int nq=4;
  if (BNEQD->QuantityFlags & QFLAG_PABS) 
   fprintf(f,"# %i absorbed power spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_PRAD) 
   fprintf(f,"# %i radiated power spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_XFORCE) 
   fprintf(f,"# %i x-force flux spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_YFORCE) 
   fprintf(f,"# %i y-force flux spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_ZFORCE) 
   fprintf(f,"# %i z-force flux spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_XTORQUE) 
   fprintf(f,"# %i x-torque flux spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_YTORQUE) 
   fprintf(f,"# %i y-torque flux spectral density\n",nq++);
  if (BNEQD->QuantityFlags & QFLAG_ZTORQUE) 
   fprintf(f,"# %i z-torque flux spectral density\n",nq++);
  fclose(f);
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("After CreateBNEQData: mem=%3.1f GB",GetMemoryUsage()/1.0e9);
  return BNEQD;

}
