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
 * buff-cas3D.h  -- header file for 3D casimir module in buff-EM suite
 *
 * homer reid     -- 2/2012
 */
#ifndef BUFFCAS3D_H
#define BUFFCAS3D_H

#include <libhrutil.h>
#include <libhmat.h>
#include <libbuff.h>

using namespace buff;

/***************************************************************/
/* brillouin zone integration schemes **************************/
/***************************************************************/
#define QUANTITY_ENERGY  1
#define QUANTITY_XFORCE  2
#define QUANTITY_YFORCE  4
#define QUANTITY_ZFORCE  8
#define QUANTITY_TORQUE1 16
#define QUANTITY_TORQUE2 32
#define QUANTITY_TORQUE3 64

#define PREAMBLE_OUT   0
#define PREAMBLE_BYXI  1

/******************************************************************/
/* BC3Data ('buff-cas3D data') is a structure that contains all   */
/* information needed to compute the contribution of a single     */
/* imaginary frequency to the Casimir quantities.                 */
/******************************************************************/
typedef struct BC3Data
 {
   SWGGeometry *G;

   int N, N1;
   HMatrix **TInvBlocks, **UBlocks, ***dUBlocks, *M, *dM;
   HVector *MInfLUDiagonal;
   int *ipiv;

   int WhichQuantities;
   int NumQuantities;

   int NumTorqueAxes;      // this number is in the range 0--3
   double TorqueAxes[9];   // [0,1,2] = x,y,z comps of 1st torque axis; [3,4,5] = 2nd axis, etc
   double GammaMatrix[27];

   GTComplex **GTCList;
   int NumTransformations;

   int NTNQ;

   int *Converged;

   double Xi;

   char *ByXiFileName;
   char *WriteCache;

   // these quantities are used to set limits for adaptive frequency integrations
   int MaxXiPoints;
   double AbsTol, RelTol;

   bool UseExistingData;

 } BC3Data;

BC3Data *CreateBC3Data(SWGGeometry *G, char *TransFile,
                       int WhichQuantities, int NumQuantities,
                       int NumTorqueAxes, double TorqueAxes[9]);

void WriteFilePreamble(FILE *f, BC3Data *BC3D, int PreambleType);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int CacheRead(const char *ByXiFileName, BC3Data *BC3D, double Xi, double *EFT);
void GetXiIntegrand(BC3Data *BC3D, double Xi, double *EFT);
void GetXiIntegral_Fixed(BC3Data *BC3D, int NumIntervals, double *EFT, double *Error);
void GetXiIntegral_Adaptive(BC3Data *BC3D, double *I, double *E);
void GetMatsubaraSum(BC3Data *BC3D, double Temperature, double *EFT, double *Error);

#endif // #define BUFFCAS3D_H
