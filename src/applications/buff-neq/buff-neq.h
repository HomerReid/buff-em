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

#define QFLAG_PABS      1
#define QFLAG_PRAD      2
#define QFLAG_XFORCE    4
#define QFLAG_YFORCE    8
#define QFLAG_ZFORCE   16
#define QFLAG_XTORQUE  32
#define QFLAG_YTORQUE  64
#define QFLAG_ZTORQUE 128

#define QINDEX_PABS     0
#define QINDEX_PRAD     1
#define QINDEX_XFORCE   2
#define QINDEX_YFORCE   3
#define QINDEX_ZFORCE   4
#define QINDEX_XTORQUE  5
#define QINDEX_YTORQUE  6
#define QINDEX_ZTORQUE  7

#define MAXQUANTITIES 8

// maximum possible number of ways to do PFT computations
#define MAXPFT 5

// methods for evaluating frequency integrals
#define QMETHOD_ADAPTIVE 1
#define QMETHOD_TRAPSIMP 2

// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in our internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

// for example: suppose in the real world we have 
// omega=3e14 rad/sec, T = 300 kelvin. then
// \hbar \omega / (kT) 
//   = (6.6e-16 ev s )(3e14 s^{-1}) / (0.026 ev) 
//   = 7.6
// whereas in this code we would have Omega==1, T=300, and hence
// Omega/(BOLTZMANNK*T) = (1/(4.36763e-4*300)) = 7.6. 

// \hbar * \omega_0^2 in units of watts
#define HBAROMEGA02 9.491145534e-06

/****************************************************************/
/* BNEQData ('buff-neq data') is a structure that contains all */
/* information needed to run computations a given frequency.    */
/****************************************************************/
typedef struct BNEQData
 {
   SWGGeometry *G;

   SVTensor **TemperatureSVTs; // one for each object in G
   double *TAvg;               // one for each object in G
   double TEnvironment;

   GTComplex **GTCList;
   int NumTransformations;

   int QuantityFlags;
   int NQ;
   int NONQ;
   int NTNONQ;

   // HMatrix structures for the BEM matrix and its subblocks
   HMatrix ***GBlocks;   // G[no][nop] = G_{no, nop} block

   // SMatrix structures for sparse matrix subblocks
   // VBlocks[no] = V matrix for object #no
   // RBlocks[no] = Rytov matrix for object #no
   SMatrix **VBlocks;
   SMatrix **RBlocks;

   HMatrix *SInverse;

   // internally-stored buffers for linear algebra operations
   HMatrix *WorkMatrix[3];

   char *FileBase;
   int NumPFTMethods;
   int PFTMethods[MAXPFT];
   char *PFTFileNames[MAXPFT];
   int DSIPoints[MAXPFT];
   PFTOptions *pftOptions;
   bool DoMomentPFT;

   bool UseExistingData;

 } BNEQData;

/****************************************************************/
/****************************************************************/
/****************************************************************/
BNEQData *CreateBNEQData(char *GeoFile, char *TransFile, int QuantityFlags,
                         char *FileBase,
                         bool DoOPFT, bool DoEMTPFT, bool DoMomentPFT,
                         int DSIPoints, double DSIRadius, char *DSIMesh,
                         int DSIPoints2);

/***************************************************************/
/* routines in GetFlux.cc **************************************/
/***************************************************************/
int GetFluxIndex(BNEQData *BNEQD, int nt, int nos, int nod, int nq);
int GetNEQIIndex(BNEQData *BNEQD, int nt, int nos, int nod, int nq);
void GetFlux(BNEQData *BNEQD, cdouble Omega, double *Flux);
void GetNEQIntegrand(BNEQData *BNEQD, cdouble Omega, double *NEQIntegrand);

/***************************************************************/
/* routines in Quadrature.cc                                   */
/***************************************************************/
void EvaluateFrequencyIntegral_Adaptive(BNEQData *BNEQD,
                                        double OmegaMin, double OmegaMax,
                                        double AbsTol, double RelTol);

void EvaluateFrequencyIntegral_TrapSimp(BNEQData *BNEQD,
                                        double OmegaMin, double OmegaMax, int NumIntervals);
#endif
