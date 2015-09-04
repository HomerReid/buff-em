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
void GetMomentPFT(SWGGeometry *G, int no, cdouble Omega,
                  HVector *JVector, HMatrix *DMatrix,
                  HMatrix *PFTMatrix,
                  double QPF[3], char *FileBase);
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
/* Note: the return value is a pointer to a preallocated buffer*/
/* within the BNEQD structure                                  */
/***************************************************************/
HMatrix *ComputeDressedRytovMatrix(BNEQData *BNEQD, int SourceObject)
{
  SWGGeometry *G        = BNEQD->G;
  HMatrix ***GBlocks    = BNEQD->GBlocks;
  SMatrix **VBlocks     = BNEQD->VBlocks;
  HMatrix *SInverse     = BNEQD->SInverse;
  SMatrix **RBlocks     = BNEQD->RBlocks;
  HMatrix **WorkMatrix  = BNEQD->WorkMatrix;

  int NO               = G->NumObjects;
  int *Offset          = G->BFIndexOffset;
  int NBF              = G->TotalBFs;

  // compute DressedRytov = W* (S^{-1}*Rytov*S^{-1}) * W^\dagger
  //                   -(S^{-1}*Rytov*S^{-1})
  //  (where W = inv( 1 + S^{-1} V S^-1 G ) )
  // in a series of steps:
  // (a) M1    <- G
  // (b) M2    <- SInverse*G
  // (c) M3    <- V
  // (d) M1    <- M3*M2 = V*SInverse*G
  // (e) M3    <- SInverse*M1 = SInverse*V*SInverse*G
  // (f) M3    += 1
  // (g) M3    <- inv(M3) = W
  // (h) M1    <- Rytov
  // (i) M2    <- SInverse*M1 = SInverse*Rytov
  // (j) M1    <- M2*SInverse = SInverse*Rytov*SInverse
  // (k) M2    <- W*M1 = M3*M1
  // (l) M3    <- M2*W'

  HMatrix *M1=WorkMatrix[0];
  HMatrix *M2=WorkMatrix[1];
  HMatrix *M3=WorkMatrix[2];

  // step (a)
  for(int no=0; no<NO; no++)
   for(int nop=no; nop<NO; nop++)
    {
      M1->InsertBlock(GBlocks[no][nop], Offset[no], Offset[nop]);
      if (nop>no)
       M1->InsertBlockTranspose(GBlocks[no][nop], Offset[nop], Offset[no]);
    };

  // step (b)
  SInverse->Multiply(M1, M2);

  // step (c)
  M3->Zero();
  for(int no=0; no<NO; no++)
   M3->AddBlock(VBlocks[no], Offset[no], Offset[no]);

  // step (d)
  M3->Multiply(M2, M1);

  // step (e)
  SInverse->Multiply(M1, M3);

  // step (f) 
  for(int nbf=0; nbf<NBF; nbf++)
   M3->AddEntry(nbf, nbf, 1.0);

  // step (g) 
  M3->LUFactorize();
  M3->LUInvert();

  // step (h)
  M1->Zero();
  M1->AddBlock(RBlocks[SourceObject], Offset[SourceObject], Offset[SourceObject]);

  // step (i)
  SInverse->Multiply(M1, M2);

  // step (j)
  M2->Multiply(SInverse, M1);

  // step (k)
  M3->Multiply(M1, M2);

  // step (l)
  M2->Multiply(M3, M1, "--transB C");

  return M1;

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
  SWGGeometry *G             = BNEQD->G;
  SVTensor **TemperatureSVTs = BNEQD->TemperatureSVTs;
  double TEnvironment        = BNEQD->TEnvironment;
  HMatrix ***GBlocks         = BNEQD->GBlocks;
  SMatrix **VBlocks          = BNEQD->VBlocks;
  SMatrix **RBlocks          = BNEQD->RBlocks;
  int NumPFTMethods          = BNEQD->NumPFTMethods;
  int *PFTMethods            = BNEQD->PFTMethods;
  char **PFTFileNames        = BNEQD->PFTFileNames;
  int *DSIPoints             = BNEQD->DSIPoints;
  PFTOptions *pftOptions     = BNEQD->pftOptions;
  char *FileBase             = BNEQD->FileBase;

  /***************************************************************/
  /* compute transformation-independent matrix blocks            */
  /***************************************************************/
  int NO=G->NumObjects;
  for(int no=0; no<NO; no++)
   { 
     Log(" Assembling V_{%i} and Rytov_{%i} ...",no,no);
     G->AssembleOverlapBlocks(no, Omega,
                              TemperatureSVTs[no], TEnvironment,
                              VBlocks[no], 0, RBlocks[no]);

     if (G->Mate[no]!=-1)
      { Log(" Block %i is identical to %i (reusing GSelf)",no,G->Mate[no]);
        continue;
      };
     Log(" Assembling G_{%i,%i}...",no,no);
     G->AssembleGBlock(no, no, Omega, GBlocks[no][no]);
   };

  static bool FirstTime=true;
  if (FirstTime)
   { FirstTime=false;
     for(int no=0; no<NO; no++) 
      { VBlocks[no]->EndAssembly();
        RBlocks[no]->EndAssembly();
      };
   };

  double PreFactor=1.0;
  for(int no=0; no<NO; no++)
   if (TemperatureSVTs[no])
    PreFactor=HBAROMEGA02;

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

     /*--------------------------------------------------------------*/
     /*- compute the requested quantities for all objects           -*/
     /*- note: nos = 'num object, source'                           -*/
     /*-       nod = 'num object, destination'                      -*/
     /*--------------------------------------------------------------*/
     static HMatrix *PFTMatrix = new HMatrix(NO, NUMPFT);
     for(int nos=0; nos<NO; nos++)
      { 
        // get the DressedRytov matrix for source object #nos
        pftOptions->RytovMatrix=ComputeDressedRytovMatrix(BNEQD, nos);

        // get the PFT for all destination objects using all
        // requested methods
        // note: the Flux vector will get filled in with the 
        // results of whichever PFT computation comes last   
        for(int nPFT=0; nPFT<NumPFTMethods; nPFT++)
         {
           pftOptions->PFTMethod  = PFTMethods[nPFT];
           pftOptions->DSIPoints  = DSIPoints[nPFT];
           G->GetPFT(0, Omega, PFTMatrix, pftOptions);
           PFTMatrix->Scale(PreFactor);

           FILE *f=fopen(PFTFileNames[nPFT],"a");
           for(int nod=0; nod<NO; nod++)
            { fprintf(f,"%s %e %i%i ",Tag,real(Omega),nos+1,nod+1);
              for(int nq=0; nq<NUMPFT; nq++) 
               fprintf(f,"%e ",PFTMatrix->GetEntryD(nod,nq));
              fprintf(f,"\n");

              for(int nq=0; nq<BNEQD->NQ; nq++)
               { int Index=GetIndex(BNEQD, nt, nos, nod, nq);
                 Flux[Index]=PFTMatrix->GetEntryD(nod,nq);
               };
            };
           fclose(f);
         };
      }; 

     /*--------------------------------------------------------------*/
     /* untransform the geometry                                     */
     /*--------------------------------------------------------------*/
     G->UnTransform();
     Log(" ...done!");

   }; // for (nt=0; nt<BNEQD->NumTransformations... )

} 
