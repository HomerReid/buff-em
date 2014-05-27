/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * TGFluctuations.h -- header file for TGFluctuations.cc
 *
 * homer reid       -- 5/2014
 */
#ifndef TGFLUCTUATIONS_H
#define TGFLUCTUATIONS_H 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libSWG.h"

using namespace scuff;
using namespace libSWG;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct TGData
 {
   int NumObjects;
   SWGVolume **Objects;
   IHAIMatProp **MPs;

   HMatrix **TMatrices;
   HMatrix ***UMatrices;

 } TGData;

#endif //#ifndef TGFLUCTUATIONS_H
