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
 * RHSVector.cc -- libSWG routines for assembling the RHS vector
 *
 * homer reid   -- 6/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libbuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

using namespace scuff;

namespace buff {

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *SWGGeometry::AllocateRHSVector()
{
  return new HVector(TotalBFs, LHM_COMPLEX);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct RHSVectorIntegrandData 
 {
   IncField *IF;
 } RHSVectorIntegrandData;

void RHSVectorIntegrand(double *x, double *b, double Divb, 
                        void *UserData, double *I)
{
  RHSVectorIntegrandData *Data = (RHSVectorIntegrandData *)UserData;

  cdouble EH[6];
  Data->IF->GetFields(x, EH);

  cdouble *zI = (cdouble *)I;
  zI[0] = b[0]*EH[0] + b[1]*EH[1] + b[2]*EH[2];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *SWGGeometry::AssembleRHSVector(cdouble Omega, IncField *IF, HVector *V)
{
  IF->SetFrequency(Omega, true);

  RHSVectorIntegrandData MyData, *Data=&MyData;
  Data->IF = IF;
 
  cdouble PreFactor = -1.0 / (II*Omega*ZVAC);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumTasks, NumThreads = GetNumThreads();
  Log("Assembling RHS vector at Omega=%g\n",z2s(Omega));
#ifndef USE_OPENMP
  NumTasks=NumThreads=1;
  Log(" no multithreading...");
#else
  NumTasks=NumThreads*100;
  Log(" OpenMP multithreading (%i threads,%i tasks)...",NumThreads,NumTasks);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int no=0; no<NumObjects; no++)
   { 
     SWGVolume *O = Objects[no];
     int Offset   = BFIndexOffset[no];
     for(int nf=0; nf<O->NumInteriorFaces; nf++)
      { 
        cdouble Entry;

        BFInt(O, nf, RHSVectorIntegrand, (void *)Data,
              2, (double *)&Entry, 0, 33, 0, 0);

        V->SetEntry( Offset + nf, PreFactor * Entry);

      }; // for(int nf=0; nf<O->NumInteriorFaces; no++, nbf++)

   }; // for(int no=0...)
  
  return V;


}

} // namespace buff
