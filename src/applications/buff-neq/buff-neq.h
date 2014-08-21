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
 * buff-neq   -- a standalone code within the buff-em suite
 *            -- for implementing the numerical TG
 *            -- approach to nonequilibrium phenomena (more 
 *            -- specifically, for computing heat radiation, 
 *            -- heat transfer, and nonequilibrium casimir forces) 
 *
 * homer reid  -- 5/2012
 */
#ifndef BUFFNEQ_H
#define BUFFNEQ_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libbuff.h>

using namespace buff;

// these are 'quantity flags' and 'quantity indices' used to 
// differentiate the various quantities that may be computed
// (power flux and i-directed momentum flux for i=x,y,z)

#define QFLAG_POWER    1
#define QFLAG_XFORCE   2
#define QFLAG_YFORCE   4
#define QFLAG_ZFORCE   8
#define QFLAG_XTORQUE 16
#define QFLAG_YTORQUE 32
#define QFLAG_ZTORQUE 64

#define QINDEX_POWER   0
#define QINDEX_XFORCE  1
#define QINDEX_YFORCE  2
#define QINDEX_ZFORCE  3
#define QINDEX_XTORQUE 4
#define QINDEX_YTORQUE 5
#define QINDEX_ZTORQUE 6

#define MAXQUANTITIES 7

#define NUMSIPFT MAXQUANTITIES

/****************************************************************/
/* BNEQData ('buff-neq data') is a structure that contains all */
/* information needed to run computations a given frequency.    */
/****************************************************************/
typedef struct BNEQData
 {
   SWGGeometry *G;
   char *WriteCache;

   int PlotFlux;

   GTComplex **GTCList;
   int NumTransformations;

   int QuantityFlags;
   int NQ;
   int NSNQ;
   int NTNSNQ;

   // HMatrix structures for the BEM matrix and its subblocks
   HMatrix ***GBlocks;     // G[ns] = G-matrix block for surface #ns
   HMatrix ***GradGBlocks; // T[ns] = G-matrix block for surface #ns

   // SMatrix structures for sparse matrix subblocks
   SMatrix **ImEps;
   SMatrix **VInv;

   // internally-stored buffers for linear algebra operations
   HMatrix **X;
   HMatrix *A, *XA;

   // Buffer[0..N] are pointers into an internally-allocated
   // chunk of memory used as workspace for computing trace formulas.
   void *Buffer[MAXQUANTITIES+1];

   char *FileBase;

 } BNEQData;

BNEQData *CreateBNEQData(char *GeoFile, char *TransFile,
                         int WhichQuantities, int PlotFlux,
                         char *FileBase, bool SymGDest);

int GetIndex(BNEQData *BNEQD, int nt, int nss, int nsd, int nq);
void GetFlux(BNEQData *BNEQD, cdouble Omega, double *Flux);

void EvaluateFrequencyIntegral(BNEQData *BNEQD,
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment,
                               double AbsTol, double RelTol,
                               double *I, double *E);

void EvaluateFrequencyIntegral2(BNEQData *BNEQD,
                                double OmegaMin, double OmegaMax,
                                double *TObjects, double TEnvironment,
                                int Intervals, double *I, double *E);

void GetSIPFTMatrices(SWGGeometry *G, int WhichObject,
                      SWGVolume *BS, double R, int NumPoints,
                      cdouble Omega, bool NeedMatrix[NUMSIPFT],
                      HMatrix *MSIPFT[NUMSIPFT]);

#endif
