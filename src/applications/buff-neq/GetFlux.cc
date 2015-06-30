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
 * GetFlux.cc  -- evaluate the matrix-trace formulas that give
 *             -- the spectral density of power/momentum flux at a 
 *             -- single frequency
 *
 * homer reid  -- 8/2014
 *
 */

#include "buff-neq.h"

#define II cdouble(0.0,1.0)

namespace buff { 

void GetDSIPFTTrace(SWGGeometry *G, cdouble Omega, HMatrix *Rytov,
                    double PFT[NUMPFT], GTransformation *GT,
                    PFTOptions *Options);
               } 

/***************************************************************/
/***************************************************************/
/***************************************************************/
GTransformation *GetFullVolumeTransformation(SWGGeometry *G,
                                             int no,
                                             bool *CreatedGT)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SWGVolume *S=G->Objects[no];
  GTransformation *GT=0;
  *CreatedGT=false;
  if ( (S->OTGT!=0) && (S->GT==0) ) 
   GT=S->OTGT;
  else if ( (S->OTGT==0) && (S->GT!=0) ) 
   GT=S->GT;
  else if ( (S->OTGT!=0) && (S->GT!=0) )
   { *CreatedGT=true;
     GT=new GTransformation(S->GT);
     GT->Transform(S->OTGT);
   };

  return GT;
}

/***************************************************************/
/* the GetFlux() routine computes a large number of quantities.*/
/* each quantity is characterized by values of four indices:   */
/*                                                             */
/*  (a) nt, the geometrical transform                          */
/*  (b) nos the source object                                  */
/*  (c) nod, the destination object                            */
/*  (d) nq, the physical quantity (i.e. power, xforce, etc.)   */
/*                                                             */
/* Given values for quantities (a)--(d), this routine computes */
/* a unique index into a big vector storing all the quantities.*/
/***************************************************************/
int GetIndex(BNEQData *BNEQD, int nt, int nos, int nod, int nq)
{
  int NO    = BNEQD->G->NumObjects;
  int NQ    = BNEQD->NQ;
  int NONQ  = NO*NQ;
  int NO2NQ = NO*NO*NQ;
  return nt*NO2NQ + nos*NONQ + nod*NQ + nq; 
}

/***************************************************************/
/* return false on failure *************************************/
/***************************************************************/
bool CacheRead(BNEQData *BNEQD, cdouble Omega, double *Flux)
{
  if (BNEQD->UseExistingData==false)
   return false;

  FILE *f=vfopen("%s.SIFlux","r",BNEQD->FileBase);
  if (!f) return false;
  Log("Attempting to cache-read flux data for Omega=%e...",real(Omega));

  int NT=BNEQD->NumTransformations;
  int NO=BNEQD->G->NumObjects;
  int NQ=BNEQD->NQ;
  int nt, nos, nod, nq;
  GTComplex **GTCList=BNEQD->GTCList;
  char *FirstTag = GTCList[0]->Tag;
  int ErrorCode, LineNum=0;
  double FileOmega, rOmega=real(Omega);

  char Line[1000];
  char *Tokens[50];
  int NumTokens, MaxTokens=50;
  bool FoundFirstTag=false;
  while( FoundFirstTag==false && fgets(Line,1000,f) )
   { 
     LineNum++;
     NumTokens=Tokenize(Line, Tokens, MaxTokens, " ");
     if (NumTokens<4) continue;
     if (strcmp(Tokens[1],FirstTag)) continue;
     sscanf(Tokens[0],"%lf",&FileOmega);
     if ( fabs(FileOmega - rOmega) > 1.0e-6*rOmega ) continue;
     FoundFirstTag=true;
   
   };
  if(FoundFirstTag==false)
   { ErrorCode=1; goto fail;}

  for(nt=0; nt<NT; nt++)
   for(nos=0; nos<NO; nos++)
    for(nod=0; nod<NO; nod++)
     { 
       if ( !(nt==0 && nos==0 && nod==0) )
        { if ( !fgets(Line,1000,f) )
           { ErrorCode=2; goto fail; }
          LineNum++;
          NumTokens=Tokenize(Line, Tokens, MaxTokens, " ");
          if ( strcmp(Tokens[1],GTCList[nt]->Tag) )
           { ErrorCode=3; goto fail; }
          sscanf(Tokens[0],"%lf",&FileOmega);
          if ( fabs(FileOmega - rOmega) > 1.0e-6*rOmega )
           { ErrorCode=4; goto fail; }
        };

       if ( NumTokens < 3+NQ ) 
        { ErrorCode=5; goto fail; }
       for(nq=0; nq<NQ; nq++)
        sscanf(Tokens[3+nq],"%le", Flux+GetIndex(BNEQD, nt, nos, nod, nq) );
     };

  // success:
   Log("...success!");
   fclose(f); 
   return true;

  fail:
   switch(ErrorCode)
    { case 1: Log("could not find first tag (fail)"); break;
      case 2: Log("line %i: unexpected end of file (fail)",LineNum); break;
      case 3: Log("line %i: wrong tag (fail)",LineNum); break;
      case 4: Log("line %i: wrong frequency (fail)",LineNum); break;
      case 5: Log("line %i: too few quantities (fail)",LineNum); break;
    };
   fclose(f); 
   return false;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFlux(BNEQData *BNEQD, cdouble Omega, double *Flux)
{
  if ( CacheRead(BNEQD, Omega, Flux) )
   return;

  Log("Computing neq quantities at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* extract fields from BNEQData structure **********************/
  /***************************************************************/
  SWGGeometry *G           = BNEQD->G;
  IHAIMatProp *Temperature = BNEQD->Temperature;
  HMatrix ***GBlocks       = BNEQD->GBlocks;
  SMatrix **VInv           = BNEQD->VInv;
  SMatrix **Sigma          = BNEQD->Sigma;
  HMatrix **WorkMatrix     = BNEQD->WorkMatrix;
  PFTOptions *pftOptions   = BNEQD->pftOptions;
  int NQ                   = BNEQD->NQ;

  /***************************************************************/
  /* compute transformation-independent matrix blocks            */
  /***************************************************************/
  int NO=G->NumObjects;
  for(int no=0; no<NO; no++)
   { 
     if (G->Mate[no]!=-1)
      { Log(" Block %i is identical to %i (reusing matrix blocks)",no,G->Mate[no]);
        continue;
      }

     Log(" Assembling G_{%i,%i}...",no,no);
     G->AssembleGBlock(no, no, Omega, GBlocks[no][no]);

     if (BNEQD->SMatricesInitialized==false)
      { VInv[no]->BeginAssembly(MAXOVERLAP);
        Sigma[no]->BeginAssembly(MAXOVERLAP);
      };

     Log(" Assembling Im V_{%i} and Sigma_{%i} ...",no,no);
     G->AssembleVInvBlock(no, Omega, Temperature, VInv[no], Sigma[no]);
     if (BNEQD->SMatricesInitialized==false)
      { VInv[no]->EndAssembly();
        Sigma[no]->EndAssembly();
      };
   };
  BNEQD->SMatricesInitialized=true;

  /***************************************************************/
  /* now loop over transformations.                              */
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  for(int nt=0; nt<BNEQD->NumTransformations; nt++)
   { 
     /*--------------------------------------------------------------*/
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     char *Tag=BNEQD->GTCList[nt]->Tag;
     G->Transform(BNEQD->GTCList[nt]);
     Log(" Computing quantities at geometrical transform %s",Tag);

     /*--------------------------------------------------------------*/
     /* assemble off-diagonal G-matrix blocks.                       */
     /* note that not all off-diagonal blocks necessarily need to    */
     /* be recomputed for all transformations; this is what the 'if' */
     /* statement here is checking for.                              */
     /*--------------------------------------------------------------*/
     for(int no=0; no<NO; no++)
      for(int nop=no+1; nop<NO; nop++)
       G->AssembleGBlock(no, nop, Omega, GBlocks[no][nop]);
       //if ( nt==0 || G->ObjectMoved[no] || G->ObjectMoved[nop] )

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     // 1. set W = G-matrix
     int *Offset = G->BFIndexOffset;
     HMatrix *W  = WorkMatrix[0];
     for(int no=0; no<NO; no++)
      for(int nop=no; nop<NO; nop++)
       { W->InsertBlock(GBlocks[no][nop], Offset[no], Offset[nop]);
         if (nop>no)
          W->InsertBlockTranspose(GBlocks[no][nop], Offset[nop], Offset[no]);
       };
     // 2. save a copy of G matrix in WGm1
     HMatrix *WGm1  = WorkMatrix[1];
     WGm1->Copy(W);
     // 3. set W = VInverse + G-matrix = VIE matrix M
     for(int no=0; no<NO; no++)
      W->AddBlock(VInv[no], Offset[no], Offset[no]);
     // 4. set WGm1 = M^{-1} * G matrix 
     W->LUFactorize();
     W->LUSolve(WGm1);
     // 5. set WGm1 = WGm1 - identity matrix
     for(int nr=0; nr<W->NR; nr++)
      WGm1->SetEntry(nr, nr, W->GetEntry(nr,nr) - 1.0 );

     /*--------------------------------------------------------------*/
     /*- compute the requested quantities for all objects           -*/
     /*- note: nos = 'num object, source'                           -*/
     /*-       nod = 'num object, destination'                      -*/
     /*--------------------------------------------------------------*/
     FILE *f=vfopen("%s.SIFlux","a",BNEQD->FileBase);
     for(int nos=0; nos<NO; nos++)
      { 
        // get the Rytov matrix for source object #nos
        HMatrix *Rytov = WorkMatrix[2];
        Rytov->Zero();
        Rytov->AddBlock(Sigma[nos], Offset[nos], Offset[nos]);
        WGm1->Multiply(Rytov, WorkMatrix[0]);
        WorkMatrix[0]->Multiply(WGm1,Rytov,"--transB C");

        // get the PFT on all destination objects #nod
        for(int nod=0; nod<NO; nod++)
         { 
           bool CreatedGT=false;
           GTransformation *FullGT
            =GetFullVolumeTransformation(G, nod, &CreatedGT);

           double PFT[MAXQUANTITIES];
           GetDSIPFTTrace(G, Omega, Rytov, PFT, FullGT, pftOptions);
           if (CreatedGT) delete FullGT;

           fprintf(f,"%s %e %i%i",Tag,real(Omega),nos+1,nod+1);
           for(int nq=0; nq<NQ; nq++)
            { fprintf(f,"%e ",PFT[nq]);
              int Index=GetIndex(BNEQD, nt, nos, nod, nq);
              Flux[Index]=PFT[nq];
            };
           fprintf(f,"\n");
         };
      };
     fclose(f);

     /*--------------------------------------------------------------*/
     /* and untransform the geometry                                 */
     /*--------------------------------------------------------------*/
     G->UnTransform();
     Log(" ...done!");

   }; // for (nt=0; nt<BNEQD->NumTransformations... )

} 
