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
                 HVector *JVector, HMatrix *Rytov,
                 HMatrix *PFTMatrix)

// PFT by displaced-surface-integral method
void GetDSIPFT(SWGGeometry *G, IncField *IF, HVector *J,
               cdouble Omega, double PFT[NUMPFT],
               GTransformation *GT, PFTOptions *Options);

void GetDSIPFTTrace(SWGGeometry *G, HMatrix *Rytov,
                    cdouble Omega, double PFT[NUMPFT],
                    GTransformation *GT, PFTOptions *Options);

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetJJ(HVector *JVector, HMatrix *Rytov, int nbfa, int nbfb)
{
  if (RytovMatrix)
   return RytovMatrix->GetEntry(nbfb, nbfa);
  else
   return conj ( JVector->GetEntry(nbfa) )
              *( JVector->GetEntry(nbfb) );
}

/***************************************************************/
/* Get the full GTransformation that takes an object/surface   */
/* from its native configuration (as described in its .msh     */
/* file) to its current configuration in a scuff-em calculation.*/
/* This is a composition of 0, 1, or 2 GTransformations        */
/* depending on whether (a) the object was DISPLACED/ROTATED   */
/* in the .scuffgeo file, and (b) the object has been          */
/* transformed since it was read in from that file.            */
/***************************************************************/
GTransformation *GetFullSurfaceTransformation(RWGGeometry *G,
                                              int ns,
                                              bool *CreatedGT)
{
  /*--------------------------------------------------------------*/
  /*- If the surface in question has been transformed since we   -*/
  /*- read it in from the meshfile (including any "one-time"     -*/
  /*- transformation specified in the .scuffgeo file) we need    -*/
  /*- to transform the cubature rule accordingly.                -*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=G->Surfaces[ns];
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
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::GetPFT(IncField *IF, HVector *JVector,
                             cdouble Omega, HMatrix *PFTMatrix,
                             PFTOptions *Options)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions DefaultOptions;
  if (Options==0)
   { Options=&DefaultOptions;
     InitPFTOptions(Options);
   };
  int PFTMethod        = Options->PFTMethod;
  HMatrix *RytovMatrix = Options->RytovMatrix; 
  HVector *RHSVector   = Options->RHSVector;

  /***************************************************************/
  /* hand off to the individual PFT algorithms to do the         */
  /* computation                                                 */
  /***************************************************************/
  if ( PFTMethod==SCUFF_PFT_OVERLAP )
   return OPFT(this, Omega, JVector, Rytov, PFTMatrix);
  else if (PFTMethod==SCUFF_PFT_EP)
   return JDEPFT(this, Omega, IF, JVector, RHSVector, 
                 Rytov, PFTMatrix);

  else if (     PFTMethod==SCUFF_PFT_DSI
             || PFTMethod==SCUFF_PFT_EPDSI
          )
   { 
     bool CreatedGT;
     GTransformation *GT
      =GetFullSurfaceTransformation(this, SurfaceIndex, &CreatedGT);

     char *DSIMesh      = Options->DSIMesh;
     double DSIRadius   = Options->DSIRadius;
     int DSIPoints      = Options->DSIPoints;
     bool DSIFarField   = Options->DSIFarField;
     double *kBloch     = Options->kBloch;
     HMatrix *RytovMatrix = Options->RytovMatrix;
     bool *NeedQuantity = Options->NeedQuantity;

     if (RytovMatrix==0)
      GetDSIPFT(this, Omega, kBloch, KN, IF, PFT,
                DSIMesh, DSIRadius, DSIPoints,
                false, DSIFarField, FluxFileName, GT);
     else 
      GetDSIPFTTrace(this, Omega, RytovMatrix,
                     PFT, NeedQuantity,
                     DSIMesh, DSIRadius, DSIPoints,
                     false, DSIFarField, FluxFileName, GT);

     if (CreatedGT) delete(GT);
   }
  else if (PFTMethod==SCUFF_PFT_EP)
   { 
     // EP force, torque
     int EPFTOrder=Options->EPFTOrder;
     double EPFTDelta=Options->EPFTDelta;
     GetEPFT(this, SurfaceIndex, Omega, KN, IF, PFT + 2,
             ByEdge, EPFTOrder, EPFTDelta);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMethod==SCUFF_PFT_EP             ||
       PFTMethod==SCUFF_PFT_EPOVERLAP      ||
       PFTMethod==SCUFF_PFT_EPDSI
     )
   {
     double Power[2];
     HMatrix *TInterior =  Options->TInterior;
     HMatrix *TExterior =  Options->TExterior;
     HMatrix *RytovMatrix = Options->RytovMatrix;
     GetEPP(this, SurfaceIndex, Omega, KN, RytovMatrix, Power,
            ByEdge, TInterior, TExterior);

     // replace absorbed and scattered power with EP calculations
     PFT[0] = Power[0];
     PFT[1] = Power[1];
   };


  /***************************************************************/
  /* produce flux plots if that was requested ********************/
  /***************************************************************/
  if (ByEdge)
   { 
     static const char *PFTNames[NUMPFT]
      ={"PAbs","PScat","FX","FY","FZ","TX","TY","TZ"};

     for(int nq=0; nq<NUMPFT; nq++)
      { char Tag[20];
        snprintf(Tag,20,"%s(%s)",PFTNames[nq],z2s(Omega));
        S->PlotScalarDensity(ByEdge[nq],true,FluxFileName,Tag);
      };

     free(ByEdge[0]);
     free(ByEdge);

   };
 
}

/***************************************************************/
/* routine for initializing a PFTOptions structure to default  */
/* values; creates and returns a new default structure if      */
/* called with Options=NULL or with no argument                */
/***************************************************************/
PFTOptions *InitPFTOptions(PFTOptions *Options)
{
  if (Options==0)
   Options = (PFTOptions *)mallocEC( sizeof(PFTOptions) );

  // general options
  Options->PFTMethod = SCUFF_PFT_DEFAULT;
  Options->FluxFileName=0;
  Options->RytovMatrix=0;
  Options->kBloch=0;

  // options affecting overlap PFT computation
  Options->RHSVector = 0;

  // options affecting DSI PFT computation
  Options->DSIMesh=0;
  Options->DSIRadius=10.0;
  Options->DSIPoints=302;
  Options->DSIFarField=false;
  for(int nq=0; nq<NUMPFT; nq++) 
   Options->NeedQuantity[nq]=true;
 
  // options affecting EP power computation
  Options->TInterior=0;
  Options->TExterior=0;

  // options affecting EP force / torque computation
  Options->EPFTOrder=1;
  Options->EPFTDelta=1.0e-5;

  Options->GetRegionPFTs=false;

  return Options;
}

} // namespace scuff
