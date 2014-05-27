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
 * GetFluxDensities.cc -- get densities of energy and momentum flux 
 *
 * homer reid          -- 5/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;
using namespace buff;

void GetFluxDensities(TGFData *Data, cdouble Omega, double *Densities)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumObjects      = Data->NumObjects;
  SWGVolume **Objects = Data->Objects;
  IHAIMatProp **MPs   = Data->MPs;

  HMatrix **TMatrices  = Data->TMatrices;
  HMatrix ***UMatrices = Data->UMatrices;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int no=0; no<NumObjects; no++)
   Objects[no]->ComputeTMatrix(Omega, MPs[no], TMatrices[no]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nt=0; nt<NumTransformations; nt++)
   {    
   };

}
