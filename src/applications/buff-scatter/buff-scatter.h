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
 * buff-scatter.h -- a standalone code within the buff-em suite
 *                -- for solving scattering problems
 *
 * homer reid     -- 6/2011--2/2012
 */
#ifndef BUFFSCATTER_H
#define BUFFSCATTER_H

#include <libhrutil.h>
#include <libhmat.h>
#include "libIncField.h"
#include "libbuff.h"

using namespace scuff;
using namespace buff;

/***************************************************************/
/* data structure containing everything needed to execute a    */
/* scattering calculation                                      */
/***************************************************************/
typedef struct BSData
 {
   SWGGeometry *G;
   HMatrix *M;
   HVector *RHS, *J;
   cdouble Omega;
   IncField *IF;
   void *opTable;
 } BSData;
 

/***************************************************************/
/* these are the 'output modules' that compute and process the */
/* scattered fields in various ways.                           */
/***************************************************************/
//void GetMoments(BSData *BSD, char *MomentFile);
void WritePFTFile(BSData *SSD, char *PFTFile, bool NeedQuantity[6]);
void ProcessEPFile(BSData *BSData, char *EPFileName);
void WriteMomentFile(BSData *BSD, char *FileName);

#endif
