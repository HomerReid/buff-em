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
 * CasimirIntegrand.cc    -- evaluate the casimir energy, force, and/or torque
 *                        -- integrand at a single Xi point or a single (Xi,kBloch)
 *                        -- point.
 *
 * homer reid  -- 10/2008 -- 6/2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <unistd.h>
#include <ctype.h>

#include "libbuff.h"
#include "buff-cas3D.h"

extern "C" {
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
}

using namespace buff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/* compute \log \det \{ M^{-1} MInfinity \}                    */
/***************************************************************/
double GetLNDetMInvMInf(BC3Data *BC3D)
{ 
  HMatrix *M              = BC3D->M;
  int N                   = BC3D->N;

  double LNDet=0.0;
  HVector *MInfLUDiagonal = BC3D->MInfLUDiagonal;

  for(int n=0; n<N; n++)
   LNDet+=log( fabs( MInfLUDiagonal->GetEntryD(n) / M->GetEntryD(n,n) ) );

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  // paraphrasing the physicists of the 1930s, 'just because
  // something is infinite doesn't mean that it's zero.' and yet...
  if (!isfinite(LNDet))
   LNDet=0.0;
  return -LNDet/(2.0*M_PI);

} 

/***************************************************************/
/* compute \trace \{ M^{-1} dMdAlpha\},                        */
/* where Alpha=x, y, z                                         */
/***************************************************************/
#if 0
double GetTraceMInvdM(BC3Data *BC3D, char XYZT)
{ 
  /***************************************************************/
  /* unpack fields from workspace structure **********************/
  /***************************************************************/
  SWGGeometry *G     = BC3D->G;
  HMatrix *M         = BC3D->M;
  HMatrix *dM        = BC3D->dM;
  HMatrix **dUBlocks = BC3D->dUBlocks;

  int Mu;
  switch(XYZT)
   {  case 'X': Mu   = 0; break;
      case 'Y': Mu   = 1; break;
      case 'Z': Mu   = 2; break;
      case '1': Mu   = 3; break; 
      case '2': Mu   = 4; break; 
      case '3': Mu   = 5; break; 
   };

  if ( '1'<=XYZT && XYZT<='3' )
   Log("  Computing torque about axis #%c...",XYZT);
  else
   Log("  Computing %cForce...",XYZT);

  /***************************************************************/
  /* stamp derivative blocks into dM matrix, compute M^{-1} dM,  */
  /* then sum the diagonals of the upper matrix block            */
  /***************************************************************/
  dM->Zero();
  for(int ns=1; ns<G->NumObjects; no++)
   dM->InsertBlockTranspose(dUBlocks[ 6*(ns-1) + Mu ], G->BFIndexOffset[ns], 0);

  M->LUSolve(dM);

  double Trace=0.0;
  for(int n=0; n<dM->NC; n++)
   Trace+=dM->GetEntryD(n,n);
  Trace*=2.0;

  // paraphrasing the physicists of the 1930s, 'just because
  // something is infinite doesn't mean that it's zero.' and yet...
  if (!isfinite(Trace))
   Trace=0.0;

  return -Trace/(2.0*M_PI);
} 
#endif


/***************************************************************/
/* stamp T and U blocks into the BEM matrix, then LU-factorize.*/
/***************************************************************/
void Factorize(BC3Data *BC3D)
{ 
  SWGGeometry *G = BC3D->G;
  HMatrix *M     = BC3D->M;

  /***************************************************************/
  /* stamp blocks into M matrix                                  */
  /***************************************************************/
  /* TInv blocks */
  for(int no=0; no<G->NumObjects; no++)
   { 
     int Offset=G->BFIndexOffset[no];
     M->InsertBlock(BC3D->TInvBlocks[no], Offset, Offset);
   };

  /* U blocks */
  for(int nb=0, no=0; no<G->NumObjects; no++)
   for(int nop=no+1; nop<G->NumObjects; nop++, nb++)
    { 
      int RowOffset=G->BFIndexOffset[no];
      int ColOffset=G->BFIndexOffset[nop];
      M->InsertBlock(BC3D->UBlocks[nb], RowOffset, ColOffset);
      M->InsertBlockTranspose(BC3D->UBlocks[nb], ColOffset, RowOffset);
    };

  /***************************************************************/
  /* LU factorize                                                */
  /***************************************************************/
  M->LUFactorize();

} 

/***************************************************************/
/* evaluate the casimir energy, force, and/or torque integrand */
/* at a single Xi point, possibly under multiple spatial       */
/* transformations.                                            */
/*                                                             */
/* the output vector 'EFT' stands for 'energy, force, and      */
/* torque.' the output quantities are packed into this vector  */
/* in the following order:                                     */
/*                                                             */
/*  energy  integrand (spatial transformation 1)               */
/*  xforce  integrand (spatial transformation 1)               */
/*  yforce  integrand (spatial transformation 1)               */
/*  zforce  integrand (spatial transformation 1)               */
/*  torque1 integrand (spatial transformation 1)               */
/*  torque2 integrand (spatial transformation 1)               */
/*  torque3 integrand (spatial transformation 1)               */
/*  energy  integrand (spatial transformation 2)               */
/* ...                                                         */
/*                                                             */
/* where only requested quantities are present, i.e. if the    */
/* user requested only energy and yforce then we would have    */
/*                                                             */
/*  EFT[0] = energy integrand (spatial transformation 1)       */
/*  EFT[1] = yforce integrand (spatial transformation 1)       */
/*  EFT[2] = energy integrand (spatial transformation 2)       */
/*  EFT[3] = yforce integrand (spatial transformation 2)       */
/*  ...                                                        */
/*  EFT[2*N-1] = yforce integrand (spatial transformation N)   */
/***************************************************************/
void GetXiIntegrand(BC3Data *BC3D, double Xi, double *EFT)
{ 
  /***************************************************************/
  /* attempt to bypass calculation by reading data from .byXi file */
  /***************************************************************/
  if ( CacheRead(BC3D->ByXiFileName, BC3D, Xi, EFT) )
   return; 

  SWGGeometry *G = BC3D->G;

  /***************************************************************/
  /* ObjectNeverMoved[ns] is initialized true and remains true   */
  /* as long as object #no is not displaced or rotated by any    */
  /* geometrical transformation, but switches to false as soon as*/
  /* that object is moved, and then remains false for the rest   */
  /* of the calculations done at this frequency                  */
  /***************************************************************/
  static bool *ObjectNeverMoved=0;
  if (ObjectNeverMoved==0)
   ObjectNeverMoved=(bool *)mallocEC(G->NumObjects*sizeof(bool));
  for(int no=0; no<G->NumObjects; no++)
   ObjectNeverMoved[no]=true;

  /***************************************************************/
  /* assemble T matrices                                         */
  /***************************************************************/
  cdouble Omega = cdouble(0.0, Xi);
  for(int no=0; no<G->NumObjects; no++)
   { 
     /* skip if this object is identical to a previous object */
     int nop = G->Mate[no];
     if ( nop != -1 )
      { Log("Object %i is identical to object %i (skipping)",no,nop);
        continue;
      };

     Log("Assembling TInv%i at Xi=%e...",no+1,Xi);
     AssembleTInvMatrix(G->Objects[no], Omega, BC3D->TInvBlocks[no]);

   }; // for(no=0; no<G->NumObjects; no++)

  /***************************************************************/
  /* if an energy calculation was requested, compute and save    */
  /* the diagonals of the LU factorization of the T blocks       */
  /* (which collectively constitute the diagonal of the LU       */
  /* factorization of the M_{\infinity} matrix).                 */
  /***************************************************************/
  if ( BC3D->WhichQuantities & QUANTITY_ENERGY )
   {
     HMatrix *M=BC3D->M;
     HVector *V=BC3D->MInfLUDiagonal;
     for(int no=0; no<G->NumObjects; no++)
      { 
        int Offset=G->BFIndexOffset[no];
        int NBF=G->Objects[no]->NumInteriorFaces;

        int nop=G->Mate[no];
        if ( nop != -1 )
         { 
           /* if this object has a mate, just copy the diagonals of the */
           /* mate's LU factor                                          */
           int MateOffset=G->BFIndexOffset[nop];
           memcpy(V->DV+Offset,V->DV+MateOffset,NBF*sizeof(double));
         }
        else
         { 
           /* LU factorize the matrix (non-PEC) case or cholesky-factorize */
           /* the negative of the matrix (all-PEC case).                   */
           /* KIND OF HACKY: we use the data buffer inside M as temporary  */
           /* storage for the content of T; this means that we have to     */
           /* call the lapack routines directly instead of using the nice  */
           /* wrappers provided by libhmat                                 */
           for(int nbf=0; nbf<NBF; nbf++)
            for(int nbfp=0; nbfp<NBF; nbfp++)
             M->DM[nbf + nbfp*NBF] = BC3D->TInvBlocks[no]->GetEntryD(nbf,nbfp);

           Log("LU-factorizing T%i at Xi=%g...",no+1,Xi);
           int info;
           dgetrf_(&NBF, &NBF, M->DM, &NBF, BC3D->ipiv, &info);
           if (info!=0)
            Log("...FAILED with info=%i (N=%i)",info,NBF);

           /* copy the LU diagonals into DRMInf */
           for(int nbf=0; nbf<NBF; nbf++)
            V->SetEntry(Offset+nbf, M->DM[nbf+nbf*NBF]);
         };
      };
   };
     
  /***************************************************************/
  /* for each line in the TransFile, apply the specified         */
  /* transformation, then calculate all quantities requested.    */
  /***************************************************************/
  int NQ = BC3D->NumQuantities;
  FILE *ByXiFile=0;
  if (BC3D->ByXiFileName)
   ByXiFile=fopen(BC3D->ByXiFileName,"a");
  for(int ntnq=0, nt=0; nt<BC3D->NumTransformations; nt++)
   { 
     char *Tag=BC3D->GTCList[nt]->Tag;

     /******************************************************************/
     /* skip if all quantities are already converged at this transform */
     /******************************************************************/
     int AllConverged=1;
     for(int nq=0; AllConverged==1 && nq<NQ ; nq++)
      if ( !BC3D->Converged[ ntnq + nq ] )
       AllConverged=0;
     if (AllConverged)
      { Log("All quantities already converged at Tag %s",Tag);

        for(int nq=0; nq<NQ; nq++)
         EFT[ntnq++]=0.0;

        if (ByXiFile)
         { fprintf(ByXiFile,"%s %.15e ",Tag,Xi);
           for(int nq=0; nq<NQ; nq++)
            fprintf(ByXiFile,"%.15e ",0.0);
           fprintf(ByXiFile,"\n");
           fflush(ByXiFile);
         };

        continue;
      };

     /******************************************************************/
     /* apply the geometrical transform                                */
     /******************************************************************/
     Log("Applying transform %s...",Tag);
     G->Transform( BC3D->GTCList[nt] );
/*
     for(int no=0; no<G->NumObjects; no++)
      if (G->ObjectMoved[no]) ObjectNeverMoved[no]=false;
*/

     /***************************************************************/
     /* assemble U_{a,b} blocks and dUdXYZT_{0,b} blocks            */
     /***************************************************************/
     for(int nb=0, no=0; no<G->NumObjects; no++)
      for(int nop=no+1; nop<G->NumObjects; nop++, nb++)
       { 
         /* if we already computed the interaction between objects ns  */
         /* and nsp once at this frequency, and if neither object has  */
         /* moved, then we do not need to recompute the interaction    */
         if ( nt>0 && ObjectNeverMoved[no] && ObjectNeverMoved[nop] )
          continue;

         Log(" Assembling U(%i,%i)",no,nop);
         AssembleGMatrix(G->Objects[no], G->Objects[nop],
                         Omega, BC3D->UBlocks[no]);

       };

     /***************************************************************/
     /* factorize the M matrix and compute Casimir quantities       */
     /***************************************************************/
     Factorize(BC3D);
     if ( BC3D->WhichQuantities & QUANTITY_ENERGY )
      EFT[ntnq++]=GetLNDetMInvMInf(BC3D);
#if 0
     if ( BC3D->WhichQuantities & QUANTITY_XFORCE )
      EFT[ntnq++]=GetTraceMInvdM(BC3D,'X');
     if ( BC3D->WhichQuantities & QUANTITY_YFORCE )
      EFT[ntnq++]=GetTraceMInvdM(BC3D,'Y');
     if ( BC3D->WhichQuantities & QUANTITY_ZFORCE )
      EFT[ntnq++]=GetTraceMInvdM(BC3D,'Z');
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE1 )
      EFT[ntnq++]=GetTraceMInvdM(BC3D,'1');
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE2 )
      EFT[ntnq++]=GetTraceMInvdM(BC3D,'2');
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE3 )
      EFT[ntnq++]=GetTraceMInvdM(BC3D,'3');
#endif

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     if (ByXiFile)
      { fprintf(ByXiFile,"%s %.15e ",Tag,Xi);
        for(int nq=0; nq<NQ; nq++)
         fprintf(ByXiFile,"%.15e ",EFT[ nt*NQ + nq ]);
        fprintf(ByXiFile,"\n");
        fflush(ByXiFile);
      };

     /******************************************************************/
     /* undo the geometrical transform                                 */
     /******************************************************************/
     G->UnTransform();

   }; // for(ntnq=nt=0; nt<BC3D->NumTransformations; nt++)
  if (ByXiFile) fclose(ByXiFile);

}
