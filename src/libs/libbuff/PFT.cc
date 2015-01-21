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
 * PFT.cc        -- libbuff class methods for computing power, force,
 *               -- and torque in classical deterministic scattering
 *               -- problems
 *
 * homer reid    -- 1/2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

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

#define NUMPFT 7

namespace buff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals(SWGVolume *OA, int nbfA, SWGVolume *OB, int nbfB,
                     cdouble Omega, cdouble IPFT[7])
{
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals(SWGVolume *OA, int nbfA, IncField *IF,
                     cdouble Omega, cdouble IPFT[7])
{
}

/***************************************************************/
/* PFT[no][nq] = nqth PFT quantity for noth object             */
/***************************************************************/
HMatrix *SWGGeometry::GetPFT(IncField *IF, HVector *JVector,
                             cdouble Omega, HMatrix *PFTMatrix)
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
  int NumThreads = GetNumThreads();
  for(int noA=0; noA<NumObjects; noA++)
   for(int noB=noA; noB<NumObjects; noB++)
    { 
      SWGObject *OA = Objects[noA];
      int OffsetA   = BFIndexOffset[noA];
      int NBFA      = OA->NumInteriorFaces;

      SWGObject *OB = Objects[noB];
      int OffsetB   = BFIndexOffset[noB];
      int NBFB      = OB->NumInteriorFaces;

      int NBFP;
      bool UseSymmetry = false;
      if (UseSymmetry)
       NBFP = OA==OB ? (NBFA*(NBFA+1)/2) : (NBFA*NBFB);
      else 
       NBFP = NBFA*NBFB;

      /*--------------------------------------------------------------*/
      /*- multithreaded loop over basis functions on OA, OB-----------*/
      /*--------------------------------------------------------------*/
      double P=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Tx=0.0, Ty=0.0, Tz=0.0; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:P, Fx, Fy, Fz, Tx, Ty, Tz)
#endif
      for(int nbfp=0; nbfp<NBFP; nbfp++)
       { 
         int nbfA, nbfB;
         if (UseSymmetry)
          { nbfA      = NBFP / NBFB;
            nbfB      = NBFP % NBFB;
          }
         else
          { nbfA      = NBFP / NBFB;
            nbfB      = NBFP % NBFB;
          };

         double IPFT[7];
         GetPFTIntegrals(OA, nbfA, OB, nbfB, Omega, IPFT);
         cdouble JJ = conj ( JVector->GetEntry(OffsetA + nbfA) )
                          *( JVector->GetEntry(OffsetB + nbfB) );

         P  += real( JJ * IPFT[0] );
         Fx += imag( JJ * IPFT[1] );
         Fy += imag( JJ * IPFT[2] );
         Fz += imag( JJ * IPFT[3] );
         Tx += imag( JJ * IPFT[4] );
         Ty += imag( JJ * IPFT[5] );
         Tz += imag( JJ * IPFT[6] );
       };

      /*--------------------------------------------------------------*/
      /*- accumulate PFT contributions for this pair of objects       */
      /*--------------------------------------------------------------*/
      double Factor = 0.5;
      if (noA==noB && nbfB > nbfA)
       Factor = 1.0;

      PFTMatrix->AddEntry(noA, 0, Factor * P  );
      if (noB!=noA)
       PFTMatrix->AddEntry(noB, 0, Factor * P );

      Factor/ = real(Omega);
      PFTMatrix->AddEntry(noA, 1, Factor * Fx );
      PFTMatrix->AddEntry(noA, 2, Factor * Fy );
      PFTMatrix->AddEntry(noA, 3, Factor * Fz );
      PFTMatrix->AddEntry(noA, 4, Factor * Tx );
      PFTMatrix->AddEntry(noA, 5, Factor * Ty );
      PFTMatrix->AddEntry(noA, 6, Factor * Tz );
      if (noB>noA)
       { PFTMatrix->AddEntry(noB, 1, -Factor * Fx );
         PFTMatrix->AddEntry(noB, 2, -Factor * Fy );
         PFTMatrix->AddEntry(noB, 3, -Factor * Fz );
         PFTMatrix->AddEntry(noB, 4, -Factor * Tx );
         PFTMatrix->AddEntry(noB, 5, -Factor * Ty );
         PFTMatrix->AddEntry(noB, 6, -Factor * Tz );
       };

    };  // for(int noA=0:NumObjects, noB=noA:NumObjects

  /***************************************************************/
  /* add incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   { 
      for(int no=0; no<NO; no++)
       { 
         SWGObject *O = Objects[no];
         int Offset   = BFIndexOffset[no];
         int NBF      = O->NumInteriorFaces;

         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         double P=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Tx=0.0, Ty=0.0, Tz=0.0; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:P, Fx, Fy, Fz, Tx, Ty, Tz)
#endif
         for(int nbf=0; nbf<O->NumInteriorFaces; nbf++)
          { 
            double IPFT[7];
            GetPFTIntegrals(O, nbf, IF, Omega, IPFT);
            cdouble J = JVector->GetEntry(Offset + nbf);
            P  += real( J*IPFT[0] );
            Fx += imag( J*IPFT[1] );
            Fy += imag( J*IPFT[2] );
            Fz += imag( J*IPFT[3] );
            Tx += imag( J*IPFT[4] );
            Ty += imag( J*IPFT[5] );
            Tz += imag( J*IPFT[6] );
          };

         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         double Factor = 0.5;
         PFTMatrix->AddEntry(no, 0, Factor * P );
         Factor/=Omega;
         PFTMatrix->AddEntry(no, 1, Factor * Fx );
         PFTMatrix->AddEntry(no, 2, Factor * Fy );
         PFTMatrix->AddEntry(no, 3, Factor * Fz );
         PFTMatrix->AddEntry(no, 4, Factor * Tx );
         PFTMatrix->AddEntry(no, 5, Factor * Ty );
         PFTMatrix->AddEntry(no, 6, Factor * Tz );
       };

   }; // if (IF)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return PFTMatrix;

}
  
} // namespace buff
