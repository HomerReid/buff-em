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

/****************************************************************/
/* BNEQData ('buff-neq data') is a structure that contains all */
/* information needed to run computations a given frequency.    */
/****************************************************************/
typedef struct BNEQData
 {
   SWGGeometry *G;

   IHAIMatProp *Temperature;

   GTComplex **GTCList;
   int NumTransformations;

   int QuantityFlags;
   int NQ;
   int NONQ;
   int NTNONQ;

   // HMatrix structures for the BEM matrix and its subblocks
   HMatrix **GBlocks;   // G[nb]     = nbth block of G-matrix
   HMatrix ***dGBlocks; // dG[nb][i] = nbth block \partial_i G-matrix

   // SMatrix structures for sparse matrix subblocks
   SMatrix **VInv;
   SMatrix **ImEps;
   bool FirstTime;

   // internally-stored buffers for linear algebra operations
   HMatrix *A, *Buf1, *Buf2;

   char *FileBase;

 } BNEQData;

#endif
