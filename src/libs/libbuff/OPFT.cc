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
 * OPFT.cc    -- libbuff class methods for computing power, force,
 *            -- and torque in classical deterministic scattering
 *            -- problems using the "overlap" formalism
 
 * homer reid -- 1/2015
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
#include "PFTOptions.h"

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
/***************************************************************/
/***************************************************************/
cdouble GetJJ(HVector *JVector, HMatrix *Rytov, int nbfa, int nbfb);

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetOPFT(SWGGeometry *G, cdouble Omega,
                 HVector *JVector, HMatrix *Rytov,
                 HMatrix *PFTMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NO = G->NumObjects;
  if (    (PFTMatrix==0)
       || (PFTMatrix->NR != NO)
       || (PFTMatrix->NC != NUMPFT)
     )
   ErrExit("invalid PFTMatrix in GetOPFT");
  PFTMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int no=0; no<NO; no++)
   { 
     SWGVolume *O      = G->Objects[no];
     int Offset        = G->BFIndexOffset[no];
     int NBF           = O->NumInteriorFaces;

     double XTorque[3] = {0.0, 0.0, 0.0};
     if (O->OTGT) O->OTGT->Apply(XTorque);
     if (O->GT)   O->GT->Apply(XTorque);

     for(int nbfA=0; nbfA<NBF; nbfA++)
      { 
        int nbfBList[MAXOVERLAP];
        double OPFTIntegrals[14][MAXOVERLAP];
        int NNZ=GetOverlapEntries(O, nbfA, OverlapIntegrand_PFT,
                                  14, (void *)XTorque,
                                  Omega, nbfBList, OPFTIntegrals);

        for(int nnz=0; nnz<NNZ; nnz++)
         { 
           int nbfB = nbfBList[nnz];
           cdouble JJ=GetJJ(JVector, Rytov, Offset+nbfA, Offset+nbfB);
           cdouble ME(OPFTIntegrals[0][nnz], OPFTIntegrals[1][nnz]);
           PFTMatrix->AddEntry(no, PFT_PABS, real(JJ*ME));

           for(int nq=PFT_XFORCE; nq<=PFT_ZTORQUE; nq++)
            { ME=cdouble( OPFTIntegrals[2*(nq-1) + 0][nnz],
                          OPFTIntegrals[2*(nq-1) + 1][nnz]
                        );
              PFTMatrix->AddEntry(no, nq, imag(JJ*ME) );
            };
         };

      }; // for(int nbfA=0; nbfA<NBF; nbfA++)

   }; // for(int no=0; no<NumObjects; no++)

  return PFTMatrix;
}
  
} // namespace buff
