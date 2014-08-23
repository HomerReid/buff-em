/*
 * buff-neq0 -- simplified version of the full buff-neq that
 *           -- computes just power, z-force, and z-torque
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "buff-neq.h"
#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
BNEQData *CreateBNEQData(SWGGeometry *G)
{
  BNEQData *BNEQD = (BNEQData *)mallocEC(sizeof(BNEQData));

  BNEQD->G = G;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NO = G->NumObjects;
  int NumBlocks = NO*(NO+1)/2;
  BNEQD->GBlocks = (HMatrix **)mallocEC(NumBlocks*sizeof(HMatrix *));
  BNEQD->dGBlocks = (HMatrix ***)mallocEC(NumBlocks*sizeof(HMatrix **));
  for(int no=0, nb=0; no<NO; no++)
   for(int nop=no; nop<NO; nop++, nb++)
    { 
      int NR=G->Objects[no]->NumInteriorFaces;
      int NC=G->Objects[nop]->NumInteriorFaces;
      BNEQD->GBlocks[nb] = new HMatrix(NR, NC, LHM_COMPLEX);
      BNEQD->dGBlocks[nb] = (HMatrix **)mallocEC(6*sizeof(HMatrix *));
      BNEQD->dGBlocks[nb][2] = new HMatrix(NR, NC, LHM_COMPLEX);
      BNEQD->dGBlocks[nb][5] = new HMatrix(NR, NC, LHM_COMPLEX);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  BNEQD->VInv  = (SMatrix **)mallocEC(NO * sizeof(SMatrix *));
  BNEQD->ImEps = (SMatrix **)mallocEC(NO * sizeof(SMatrix *));
  for(int no=0; no<NO; no++)
   {
     int NBF=G->Objects[no]->NumInteriorFaces;

     BNEQD->VInv[no] = new SMatrix(NBF, NBF, LHM_COMPLEX);
     BNEQD->VInv[no]->BeginAssembly(7*NBF);

     BNEQD->ImEps[no] = new SMatrix(NBF, NBF, LHM_REAL);
     BNEQD->ImEps[no]->BeginAssembly(7*NBF);
   };
  BNEQD->FirstTime=true;

  /***************************************************************/
  /* G, W, A are all needed simultaneously in ComputeAMatrix()   */
  /* X, A, and XA are all needed simultaneously in GetFlux()     */
  /*  --> need A matrix plus two scratch buffers to store W,G    */
  /*  --> or X,XA.                                               */
  /***************************************************************/
  int N = G->TotalBFs;
  BNEQD->A    = new HMatrix(N, N, LHM_COMPLEX);
  BNEQD->Buf1 = new HMatrix(N, N, LHM_COMPLEX);
  BNEQD->Buf2 = new HMatrix(N, N, LHM_COMPLEX);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFlux(BNEQData *BNEQD, int Quantity, cdouble *Fluxes)
                                            //double **Fluxes)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G  = BNEQD->G;

  HMatrix *A          = BNEQD->A;
  HMatrix *X          = BNEQD->Buf1;
  HMatrix *XA         = BNEQD->Buf2;
  HMatrix **GBlocks   = BNEQD->GBlocks;
  HMatrix ***dGBlocks = BNEQD->dGBlocks;

  int NO = G->NumObjects;

  /***************************************************************/
  /* assemble X matrix and compute X*A matrix                    */
  /***************************************************************/
  for(int no=0, nb=0; no<NO; no++)
   for(int nop=no+1; nop<NO; nop++, nb++)
    { 
      HMatrix *XBlock 
       = (Quantity==QINDEX_POWER) ? GBlocks[nb] : dGBlocks[nb][Quantity-1];

      int RowOffset = G->BFIndexOffset[no];
      int ColOffset = G->BFIndexOffset[nop];
      X->InsertBlock(XBlock, RowOffset, ColOffset);
      if(nop>no) 
       X->InsertBlockTranspose(XBlock, ColOffset, RowOffset);

    };
  X->Multiply(A, XA);

  /***************************************************************/
  /* loop over source and destination objects                    */
  /***************************************************************/
  for(int ns=0; ns<NO; ns++)
   for(int nd=0; nd<NO; nd++)
    { 
      SMatrix *ImEps = BNEQD->ImEps[ns];

      int sLength = G->Objects[ns]->NumInteriorFaces;
      int sOffset = G->BFIndexOffset[ns];

      int dLength = G->Objects[nd]->NumInteriorFaces;
      int dOffset = G->BFIndexOffset[nd];

      // loop over nonzero entries in the ImEps matrix
      // to compute the sum
      //  Flux = \sum_{ijk} (XA)_{ij} ImEps_{jk} A^\dagger_{ki}
      //       = \sum_{jk} ImEps_{jk} * \sum_{i} (XA)_{ij} A*_{ik}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
      //double Flux = 0.0;
      cdouble Flux = 0.0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
      for(int j=0; j<sLength; j++)
       { 
         int *nc;
         double *Entries;
         int NNZ = ImEps->GetRow(j, &nc, (void **)&Entries);

         for(int nnz=0; nnz<NNZ; nnz++)
          {
            int k = nc[nnz];
            cdouble Sum=0.0;
            for(int i=0; i<dLength; i++)
             Sum += XA->GetEntry(dOffset+i, sOffset+j) * conj(A->GetEntry(dOffset+i, sOffset+k));

            Flux += Entries[nnz] * Sum;
          };
       };

      Fluxes[ NO*ns + nd ] = Flux;
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeAMatrix(BNEQData *BNEQD)
{
  SWGGeometry *G    = BNEQD->G;
  HMatrix *A        = BNEQD->A;
  HMatrix *GG       = BNEQD->Buf1;
  HMatrix *M        = BNEQD->Buf2;
  HMatrix **GBlocks = BNEQD->GBlocks;
  SMatrix **VInv    = BNEQD->VInv;

  /*******************************************************************/
  /* assemble G matrix                                               */
  /*******************************************************************/ 
  int NO = G->NumObjects;
  for(int no=0, nb=0; no<NO; no++)
   { 
     int RowOffset = G->BFIndexOffset[no];
     GG->InsertBlock(GBlocks[nb], RowOffset, RowOffset);

     for(int nop=no+1; nop<NO; nop++, nb++)
      {    
        int ColOffset = G->BFIndexOffset[no];
        GG->InsertBlock(GBlocks[nb], RowOffset, ColOffset);
        GG->InsertBlockTranspose(GBlocks[nb], ColOffset, RowOffset);
      };
   };

  /*******************************************************************/
  /* assemble M=VInv + G matrix, then solve the equation M*(WG) = G  */
  /* to compute WG [ where W=M^{-1} ] and then A=1-WG                */
  /*******************************************************************/
  M->Copy(GG);
  for(int no=0, nb=0; no<NO; no++)
   { 
     int RowOffset = G->BFIndexOffset[no];
     M->AddBlock(VInv[no], RowOffset, RowOffset);
   };
  M->LUFactorize();

  A->Copy(GG);
  M->LUSolve(A); // A <- W*G 

  for(int nr=0; nr<A->NR; nr++)
   for(int nc=0; nc<A->NC; nc++)
    A->SetEntry(nr, nc, (nr==nc ? 1.0 : 0.0) - A->GetEntry(nr,nc) );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *OmegaFile=0;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
/**/     
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  0,  "list of (angular) frequencies"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*******************************************************************/
  /* read in the geometry ********************************************/
  /*******************************************************************/
  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  SWGGeometry *G = new SWGGeometry(GeoFile);

  int NO = G->NumObjects;

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaList=new HVector(OmegaFile,LHM_TEXT);
  if (OmegaList->ErrMsg)
   ErrExit(OmegaList->ErrMsg);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  BNEQData *BNEQD=CreateBNEQData(G);

  HMatrix **GBlocks = BNEQD->GBlocks;
  HMatrix ***dGBlocks = BNEQD->dGBlocks;
  SMatrix **VInv  = BNEQD->VInv;
  SMatrix **ImEps = BNEQD->ImEps;

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  cdouble PowerFlux[NO*NO];
  cdouble ZForceFlux[NO*NO];
  cdouble ZTorqueFlux[NO*NO];
  FILE *f=vfopen("%s.flux","w",GetFileBase(G->GeoFileName));
  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
   { 
     cdouble Omega = OmegaList->GetEntry(nFreq);

     /*******************************************************************/
     /* assemble transformation-independent matrix blocks ***************/
     /*******************************************************************/
     for(int no=0; no<NO; no++)
      { int nb = no*NO - no*(no-1)/2;
        G->AssembleGBlock(no, no, Omega, GBlocks[nb], dGBlocks[nb]);
        G->AssembleVInvBlock(no, Omega, VInv[no], ImEps[no]);
      };

     if (BNEQD->FirstTime)
      { for(int no=0; no<NO; no++)
         { VInv[no]->EndAssembly();
           ImEps[no]->EndAssembly();
         };
        BNEQD->FirstTime=false;
      };
      
     /*******************************************************************/
     /* assemble transformation-dependent matrix blocks *****************/
     /*******************************************************************/
     for(int no=0, nb=0; no<NO; no++)
      for(int nop=no+1; nop<NO; nop++, nb++)
       { 
         // insert logic to avoid recomputation as necessary 
         G->AssembleGBlock(no, nop, Omega, GBlocks[nb], dGBlocks[nb]);
       };

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     // compute A = (1-WG)
     ComputeAMatrix(BNEQD);

     // evaluate traces of all requested quantities
     GetFlux(BNEQD, QINDEX_POWER,   PowerFlux);
     GetFlux(BNEQD, QINDEX_ZFORCE,  ZForceFlux);
     GetFlux(BNEQD, QINDEX_ZFORCE,  ZTorqueFlux);

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     fprintf(f,"%e ",real(Omega));
     for(int no=0; no<NO; no++)
      for(int nop=no+1; nop<NO; nop++)
       { 
         fprintf(f,"%i%i ",no,nop);
         fprintf(f,"%s ",CD2S(PowerFlux[ NO*no + nop]));
         fprintf(f,"%s ",CD2S(ZForceFlux[ NO*no + nop]));
         fprintf(f,"%s ",CD2S(ZTorqueFlux[ NO*no + nop]));
       };
     fprintf(f,"\n");
     fflush(f);

   }; // for(int nFreq=0; nFreq<OmegaList->NumFreqs; nFreq++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fclose(f);
  printf("Thank you for your support.\n");

}

#if 0
0 1 2 3
  4 5 6                       
    7 8 
      9 

no | nb
-------
 0 |  0
 1 |  NO
 2 |  2*NO-1
 3 |  3*NO-3
 4 |  4*NO-6
no | no*NO - no*(no-1)/2;
#endif 
