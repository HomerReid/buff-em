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
 * GetPFT.cc -- libbuff class method for computing power, force,
 *           -- and torque in classical deterministic scattering
 *           -- problems. GetPFT() is actually just a switchboard
 *           -- routine that selects from among the three
 *           -- available methods for computing PFTs.
 *
 *           -- Note: for convenience, BUFF-EM uses the same
 *           -- PFTOptions structure that is defined in SCUFF-EM,
 *           -- even though some fields in that structure have
 *           -- different interpretations or are not used 
 *           -- by BUFF-EM.
 *
 * homer reid -- 6/2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fenv.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0,1)

using namespace scuff;

namespace buff {

/***************************************************************/
/* function prototypes for the various PFT algorithms.         */
/*                                                             */
/* Note: In an effort to improve modularity and readability,   */
/*       the specific PFT algorithms are implemented as        */
/*       standalone (non-class-method) functions.              */
/*       Only the master GetPFT() routine is a class method in */
/*       SWGGeometry; it is just a switchboard routine that    */
/*       hands off to the various non-class-method functions   */
/*       to do the computation.                                */
/***************************************************************/

// PFT by overlap method
HMatrix *GetOPFT(SWGGeometry *G, cdouble Omega,
                 HVector *JVector, HMatrix *DRMatrix,
                 HMatrix *PFTMatrix);

// PFT by energy/momentum transfer method
HMatrix *GetEMTPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HVector *RHSVector,
                   HMatrix *DRMatrix, HMatrix *PFTMatrix);

// PFT from multipole moments
void GetMomentPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                  HVector *JVector, HMatrix *DRMatrix,
                  HMatrix *PFTMatrix, bool KeepQpTerm=false,
                  char *FileBase=0);

// PFT by displaced-surface-integral method
void GetDSIPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
               HVector *JVector, double PFT[NUMPFT],
               GTransformation *GT1, GTransformation *GT2,
               PFTOptions *Options);

void GetDSIPFTTrace(SWGGeometry *G, cdouble Omega, HMatrix *DRMatrix,
                    double PFT[NUMPFT],
                    GTransformation *GT1, GTransformation *GT2,
                    PFTOptions *Options);

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetJJ(HVector *JVector, HMatrix *DRMatrix, int nbfa, int nbfb)
{
  if (DRMatrix)
   return DRMatrix->GetEntry(nbfb, nbfa);
  else
   return conj ( JVector->GetEntry(nbfa) )
              *( JVector->GetEntry(nbfb) );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
SWGVolume *ResolveNBF(SWGGeometry *G, int nbf, int *pno, int *pnf)
{
  int no=0, NOm1=G->NumObjects - 1;
  while( (no < NOm1) && (nbf >= G->BFIndexOffset[no+1]) )
   no++;
  if (pno) *pno = no;
  if (pnf) *pnf = nbf - G->BFIndexOffset[no];
  return G->Objects[no];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::GetPFT(HVector *JVector,
                             cdouble Omega,
                             HMatrix *PFTMatrix,
                             PFTOptions *Options)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMatrix!=0 && (PFTMatrix->NR!=NumObjects || PFTMatrix->NC!=NUMPFT) )
   { delete PFTMatrix;
     PFTMatrix=0;
   };
  if (PFTMatrix==0)
   PFTMatrix= new HMatrix(NumObjects, NUMPFT);

  PFTMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions DefaultOptions;
  if (Options==0)
   { Options=&DefaultOptions;
     InitPFTOptions(Options);
   };
  int PFTMethod      = Options->PFTMethod;
  HMatrix *DRMatrix  = Options->DRMatrix;
  HVector *RHSVector = Options->RHSVector;
  IncField *IF       = Options->IF;

  /***************************************************************/
  /* hand off to the individual PFT algorithms to do the         */
  /* computation                                                 */
  /***************************************************************/
  if ( PFTMethod==SCUFF_PFT_OVERLAP )
   GetOPFT(this, Omega, JVector, DRMatrix, PFTMatrix);
  else if ( PFTMethod==SCUFF_PFT_EMT )
   GetEMTPFT(this, Omega, IF, JVector, RHSVector, DRMatrix, PFTMatrix);
  else if ( PFTMethod==SCUFF_PFT_MOMENTS )
   GetMomentPFT(this, Omega, IF, JVector, DRMatrix, PFTMatrix);
  else // ( PFTMethod==SCUFF_PFT_DSI )
   for(int no=0; no<NumObjects; no++)
    { 
      GTransformation *GT1=Objects[no]->OTGT;
      GTransformation *GT2=Objects[no]->GT;
      double PFT[NUMPFT];
      if (DRMatrix==0)
       GetDSIPFT(this, Omega, IF, JVector, PFT, GT1, GT2, Options);
      else
       GetDSIPFTTrace(this, Omega, DRMatrix, PFT, GT1, GT2, Options);

      PFTMatrix->SetEntriesD(no, ":", PFT);
    };

  return PFTMatrix;
 
}

/***************************************************************/
/* routine for initializing a PFTOptions structure to default  */
/* values; creates and returns a new default structure if      */
/* called with Options=NULL or with no argument                */
/***************************************************************/
PFTOptions *BUFF_InitPFTOptions(PFTOptions *Options)
{
  if (Options==0)
   Options = (PFTOptions *)mallocEC( sizeof(PFTOptions) );

  // set SCUFF-EM defaults
  scuff::InitPFTOptions(Options);

  // modify as appropriate for BUFF-EM defaults
  Options->PFTMethod = SCUFF_PFT_EMT;

  // options affecting DSI PFT computation
  Options->DSIMesh=0;
  Options->DSIRadius=5.0;
  Options->DSIPoints=110;

  return Options;
}

} // namespace scuff
