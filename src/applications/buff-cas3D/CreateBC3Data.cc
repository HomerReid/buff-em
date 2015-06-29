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
 * CreateBC3Data.cc --a utility function to initialize a
 *                 -- buff-cas3D data structure for a given
 *                 -- run of the code
 *
 * homer reid      -- 5/2014
 *
 */

#include <time.h>
#include <sys/time.h>

#include <libhrutil.h>

#include "buff-cas3D.h"


using namespace scuff;
using namespace buff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
BC3Data *CreateBC3Data(SWGGeometry *G, char *TransFile,
                       int WhichQuantities, int NumQuantities,
                       int NumTorqueAxes, double TorqueAxes[9])
{
  BC3Data *BC3D=(BC3Data *)mallocEC(sizeof(*BC3D));
  BC3D->G = G;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  BC3D->WhichQuantities=WhichQuantities;
  BC3D->NumQuantities=NumQuantities;
  BC3D->NumTorqueAxes=NumTorqueAxes;
  if (BC3D->NumTorqueAxes)
   ErrExit("torque is not currently supported");
/*
  if (BC3D->NumTorqueAxes)
   { memcpy(BC3D->TorqueAxes, TorqueAxes, 3*NumTorqueAxes*sizeof(cdouble));
     for(int nta=0; nta<NumTorqueAxes; nta++)
      CreateGammaMatrix(TorqueAxes+3*nta, BC3D->GammaMatrix+9*nta);
   };
*/

  /*--------------------------------------------------------------*/
  /*- read the transformation file if one was specified and check */
  /*- that it plays well with the specified geometry file.        */
  /*- note if TransFile==0 then this code snippet still works; in */
  /*- this case the list of GTComplices is initialized to contain */
  /*- a single empty GTComplex and the check automatically passes.*/
  /*--------------------------------------------------------------*/
  BC3D->GTCList=ReadTransFile(TransFile, &(BC3D->NumTransformations));
  char *ErrMsg=G->CheckGTCList(BC3D->GTCList, BC3D->NumTransformations);
  if (ErrMsg)
   ErrExit("file %s: %s",TransFile,ErrMsg);

  BC3D->NTNQ = BC3D->NumTransformations * BC3D->NumQuantities;

  BC3D->Converged = (int *)mallocEC( (BC3D->NTNQ) * sizeof(int) );

  /*--------------------------------------------------------------*/
  /*- allocate arrays of matrix subblocks that allow us to reuse  */
  /*- chunks of the BEM matrices for multiple geometrical         */
  /*- transformations.                                            */
  /*-                                                             */
  /*- TInvBlocks[no]    = (no,no) (diagonal) block                */
  /*- UBlocks[nb]       = nbth above-diagonal block, where        */
  /*-                      nb==0   <--> (0,1) block               */
  /*-                      nb==1   <--> (0,2) block               */
  /*-                      ...     <--> ...   block               */
  /*-                      nb==N   <--> (0,N) block               */
  /*-                      nb==N+1 <--> (1,2) block               */
  /*-                     etc.                                    */
  /*-                                                             */
  /*- dUBlocks[nop][Mu] = Mu derivative of (0,nop) block          */
  /*-                     where Mu=0,1,2, for x,y,z-displacement  */
  /*-                           Mu=3,4,5, for axis 1,2,3 rotation */
  /*--------------------------------------------------------------*/
  int NO=G->NumObjects;
  BC3D->TInvBlocks  = (HMatrix **)mallocEC(NO*sizeof(HMatrix *));
  for(int no=0; no<G->NumObjects; no++)
   { int nop = G->Mate[no];
     if ( nop!=-1 )
      BC3D->TInvBlocks[no] = BC3D->TInvBlocks[nop];
     else
      { int NBF=G->Objects[no]->NumInteriorFaces;
        BC3D->TInvBlocks[no]
         = new HMatrix(NBF, NBF, LHM_REAL, LHM_SYMMETRIC);
      };
   };

  // BC3D->UBlocks[0]    = 0,1    block
  // BC3D->UBlocks[1]    = 0,2    block
  //             ...    = ...
  // BC3D->UBlocks[NO-1] = 0,NO   block
  // BC3D->UBlocks[NO]   = 1,2    block
  // etc.                                         
  int NumBlocks   = NO*(NO-1)/2; // number of above-diagonal blocks 
  BC3D->UBlocks   = (HMatrix **)mallocEC(   NumBlocks * sizeof(HMatrix *));

  // BC3D->dUBlocks[nop][Mu] = Mu derivative of (0,nop) block  */
  /*- where Mu=0,1,2, for x,y,z-displacement                   */
  /*-       Mu=3,4,5, for axis 1,2,3 rotation                  */
  BC3D->dUBlocks
   = (HMatrix ***) mallocEC( G->NumObjects * sizeof(HMatrix **) );
  BC3D->dUBlocks[0] = 0;
  for(int nop=1; nop<G->NumObjects; nop++)
   BC3D->dUBlocks[0] = (HMatrix **)mallocEC( 6*sizeof(HMatrix *) );

  for(int nb=0, no=0; no<NO; no++)
   for(int nop=no+1; nop<NO; nop++, nb++)
    { 
      int NBF=G->Objects[no]->NumInteriorFaces;
      int NBFp=G->Objects[nop]->NumInteriorFaces;

      BC3D->UBlocks[nb] = new HMatrix(NBF, NBFp);

      if ( no==0 )
       { if ( WhichQuantities & QUANTITY_XFORCE )
          BC3D->dUBlocks[nop][0] = new HMatrix(NBF, NBFp);
         if ( WhichQuantities & QUANTITY_YFORCE )
          BC3D->dUBlocks[nop][1] = new HMatrix(NBF, NBFp);
         if ( WhichQuantities & QUANTITY_ZFORCE )
          BC3D->dUBlocks[nop][2] = new HMatrix(NBF, NBFp);
         if (WhichQuantities & QUANTITY_TORQUE1)
          BC3D->dUBlocks[nop][3] = new HMatrix(NBF, NBFp);
         if (WhichQuantities & QUANTITY_TORQUE2)
          BC3D->dUBlocks[nop][4] = new HMatrix(NBF, NBFp);
         if (WhichQuantities & QUANTITY_TORQUE3)
          BC3D->dUBlocks[nop][5] = new HMatrix(NBF, NBFp);
       };

    };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N   = BC3D->N  = BC3D->G->TotalBFs;
  int N1  = BC3D->N1 = BC3D->G->Objects[0]->NumInteriorFaces;
  BC3D->M            = new HMatrix(N,  N);
  BC3D->dM           = new HMatrix(N,  N1);

  if (WhichQuantities & QUANTITY_ENERGY)
   { BC3D->MInfLUDiagonal = new HVector(G->TotalBFs);
     BC3D->ipiv = (int *)mallocEC(N*sizeof(int));
   }
  else
   { BC3D->MInfLUDiagonal=0;
     BC3D->ipiv=0;
   };

  return BC3D;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFilePreamble(FILE *f, BC3Data *BC3D, int PreambleType)
{
  if (f==0) return;

  char DateStr[40];
  time_t MyTime = time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(DateStr,30,"%D::%T",MyTm);
  fprintf(f,"# buff-cas3D run on %s at %s",GetHostName(),DateStr);
  fprintf(f,"# data file columns: \n");

  int nc=1;

  fprintf(f,"#%i: transform tag\n",nc++);

  if ( PreambleType==PREAMBLE_BYXI )
   fprintf(f,"#%i: imaginary angular frequency\n",nc++);

  if (PreambleType == PREAMBLE_OUT) 
   {
     if ( BC3D->WhichQuantities & QUANTITY_ENERGY )
      { fprintf(f,"#%i: energy \n",nc++);
        fprintf(f,"#%i: energy error \n",nc++);
      };
     if ( BC3D->WhichQuantities & QUANTITY_XFORCE )
      { fprintf(f,"#%i: x-force \n",nc++);
        fprintf(f,"#%i: x-force error \n",nc++);
      };
     if ( BC3D->WhichQuantities & QUANTITY_YFORCE )
      { fprintf(f,"#%i: y-force \n",nc++);
        fprintf(f,"#%i: y-force error \n",nc++);
      };
     if ( BC3D->WhichQuantities & QUANTITY_ZFORCE )
      { fprintf(f,"#%i: z-force \n",nc++);
        fprintf(f,"#%i: z-force error \n",nc++);
      };
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE1 )
      { fprintf(f,"#%i: 1-torque \n",nc++);
        fprintf(f,"#%i: 1-torque error \n",nc++);
      };
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE2 )
      { fprintf(f,"#%i: 2-torque \n",nc++);
        fprintf(f,"#%i: 2-torque error \n",nc++);
      };
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE3 )
      { fprintf(f,"#%i: 3-torque \n",nc++);
        fprintf(f,"#%i: 3-torque error \n",nc++);
      };
   }
  else
   {
     if ( BC3D->WhichQuantities & QUANTITY_ENERGY )
      fprintf(f,"#%i: energy integrand \n",nc++);
     if ( BC3D->WhichQuantities & QUANTITY_XFORCE )
      fprintf(f,"#%i: x-force integrand \n",nc++);
     if ( BC3D->WhichQuantities & QUANTITY_YFORCE )
      fprintf(f,"#%i: y-force integrand \n",nc++);
     if ( BC3D->WhichQuantities & QUANTITY_ZFORCE )
      fprintf(f,"#%i: z-force integrand \n",nc++);
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE1 )
      fprintf(f,"#%i: 1-torque integrand \n",nc++);
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE2 )
      fprintf(f,"#%i: 2-torque integrand \n",nc++);
     if ( BC3D->WhichQuantities & QUANTITY_TORQUE3 )
      fprintf(f,"#%i: 3-torque integrand \n",nc++);
   };

}
