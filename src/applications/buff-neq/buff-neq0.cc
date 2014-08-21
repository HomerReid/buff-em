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
void GetFlux(BNEQData *BNEQD, HMatrix *X, double **Fluxes)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G  = BNEQD->G;
  HMatrix *A      = BNEQD->A;
  HMatrix *XA     = BNEQD->XA;

  X->Multiply(A, XA);

  /***************************************************************/
  /* loop over source and destination objects                    */
  /***************************************************************/
  int NO = G->NumObjects;
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
      double Flux = 0.0;
      for(int j=0; j<sLength; j++)
       { 
         int *nc;
         double *Entries;
         int NNZ = Sigma->GetRow(nr, &nc, (void **)&Entries);

         for(int nnz=0; nnz<NNZ; nnz++)
          { 
            int k = nc[nnz];
            double Sum=0.0;
            for(int i=0; i<dLength; i++)
             Sum+=XA->GetEntry(dOffset+i, sOffset+j) * conj(A->GetEntry(dOffset+i, sOffset+k);

            Flux += Entries[nnz] *Sum;
          };
       };

      Fluxes[ns][nd] = Flux;
    };
}

void ComputeAMatrix(
{

  /*******************************************************************/
  /* assemble G matrix                                               */
  /*******************************************************************/
  for(int no=0, nb=0; no<NO; no++)
   { 
     int RowOffset = G->BFIndexOffset[no];
     GMatrix->InsertBlock(GBlocks[nb], RowOffset, RowOffset);

     for(int nop=no+1; nop<NO; nop++, nb++)
      {    
        int ColOffset = G->BFIndexOffset[no];
        GMatrix->InsertBlock(GBlocks[nb], RowOffset, ColOffset);
        GMatrix->InsertBlockTranspose(GBlocks[nb], ColOffset, RowOffset);
      };
   };

  /*******************************************************************/
  /* assemble M=VInv + G matrix, then invert it to get W matrix      */
  /*******************************************************************/
  W->Copy(GMatrix);
  for(int no=0, nb=0; no<NO; no++)
   { 
     int RowOffset = Geom->BFIndexOffset[no];
     W->AddBlock(VInv[no], RowOffset, RowOffset);
   };
  W->LUFactorize();
  W->LUInvert();

  /*******************************************************************/
  /* compute 1 - W*G                                                 */
  /*******************************************************************/
  W->Multiply(G, A ); // A <- W*G
  for(int nr=0; nr<A->NR; nr++)
   for(int nc=0; nc<A->NC; nc++)
    A->SetEntry(nr, nc, (nr==nc ? 1.0 : 0.0) - A->GetEntry(nr,nc));
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
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
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
  int N  = G->TotalBFs;
  int N1 = G->Objects[0]->NumInteriorFaces;
  int N2 = G->Objects[1]->NumInteriorFaces;

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  if (OmegaFile==0)
   OSUsage(argv[0], OSArray, "--OmegaFile option is mandatory");

  OmegaList=new HVector(OmegaFile,LHM_TEXT);
  if (OmegaList->ErrMsg)
   ErrExit(OmegaPoints->ErrMsg);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix **GBlocks, ***GradGBlocks;
  int NumBlocks = NO*(NO+1)/2;
  GBlocks = (HMatrix **)malloc(NumBlocks*sizeof(HMatrix *));
  dGBlocks = (HMatrix ***)malloc(NumBlocks*sizeof(HMatrix **));
  for(int no=0, nb=0; no<NO; no++)
   for(int nop=no; nop<NO; nop++, nb++)
    { 
      int NR=Geom->Objects[no]->NumInteriorFaces;
      int NC=Geom->Objects[nop]->NumInteriorFaces;
      GBlocks[nb] = new HMatrix(NR, NC, LHM_COMPLEX);
      dGBlocks[nb] = (HMatrix **)malloc(6*sizeof(HMatrix *));
      memset(dGBlocks[nb], 0, 6*sizeof(HMatrix *));
      dGBlocks[nb][2] = new HMatrix(NR, NC, LHM_COMPLEX);
      dGBlocks[nb][5] = new HMatrix(NR, NC, LHM_COMPLEX);
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SMatrix **VInv, **ImEps;
  VInv  = (SMatrix **)malloc(NO * sizeof(SMatrix *));
  ImEps = (SMatrix **)malloc(NO * sizeof(SMatrix *));
  for(int no=0; no<NO; no++)
   { int N=Geom->Objects[no]->NumInteriorFaces;
     VInv[no] = new SMatrix(N, N, LHM_COMPLEX);
     ImEps[no] = new SMatrix(N, N, LHM_REAL);
     VInv[no]->BeginAssembly(7*N);
     ImEps[no]->BeginAssembly(7*N);
   };

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  HMatrix *G, *W, *A;
  W    = new HMatrix(N, N, LHM_COMPLEX);
  IMWG = new HMatrix(N, N, LHM_COMPLEX);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  static bool FirstTime=true;
  for(int nFreq=0; nFreq<OmegaList->NumFreqs; nFreq++)
   { 
     cdouble Omega = OmegaList->GetEntry(nFreq);

     /*******************************************************************/
     /* assemble transformation-independent matrix blocks ***************/
     /*******************************************************************/
     for(int no=0; no<NO; no++)
      { int nb = no*NO - no*(no-1)/2;
        G->AssembleGBlock(no, no, Omega, GBlocks[nb], GradGBlocks[nb]);
        G->AssembleVInvBlock(no, Omega, VInv[no], ImEps[no]);

        if (FirstTime)
         { VInv[no]->EndAssembly();
           ImEps[no]->EndAssembly();
         };

      };
     FirstTime=false;
      
     /*******************************************************************/
     /* assemble transformation-dependent matrix blocks *****************/
     /*******************************************************************/
     for(int no=0, nb=0; no<NO; no++)
      for(int nop=no+1; nop<NO; nop++, nb++)
       { 
         // insert logic to avoid recomputation as necessary 
         G->AssembleGBlock(no, nop, Omega, GBlocks[nb], GradGBlocks[nb]);

         // compute A = (1-WG)
         ComputeA(BNEQD);

         // evaluate traces of all requested quantities
         GetFlux(BNEQD, QUANTITY_POWER,   PowerFlux);
         GetFlux(BNEQD, QUANTITY_ZFORCE,  ZForceFlux);
         GetFlux(BNEQD, QUANTITY_ZTORQUE, ZTorqueFlux);

       };
   { 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
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
