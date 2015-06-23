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

 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/***************************************************************/
/* DSIPFT.cc  -- BUFF-EM code for computing power, force, and  */
/*            -- torque (PFT) using the 'displaced surface     */
/*            -- integral' (DSI) method                        */
/*                                                             */
/* Homer Reid -- 6/2015                                        */
/***************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libTriInt.h>
#include <libscuff.h>
#include <PFTOptions.h>
#include <GTransformation.h>
#include <libbuff.h>

#ifdef USE_OPENMP
 #include <omp.h>
#endif
 
#define II cdouble(0.0,1.0) 

namespace scuff {

HMatrix *GetSCRMatrix(char *BSMesh, double R, int NumPoints, 
                      bool UseCCQ, GTransformation *GT);

void GetNMatrices(double nHat[3], double X[3], double XTorque[3],
                  double NMatrix[NUMPFT][3][3], 
                  bool *NeedQuantity=0);

double HVMVP(cdouble V1[3], double M[3][3], cdouble V2[3]);

                }

namespace buff{

/***************************************************************/
/* Get power, force, and torque by the displaced               */
/* surface-integral method.                                    */
/***************************************************************/
void GetDSIPFT(SWGGeometry *G, IncField *IF, HVector *J,
               cdouble Omega, double PFT[NUMPFT],
               char *DSIMesh, double DSIRadius, int DSIPoints,
               GTransformation *GT)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (DSIMesh)
   Log("Computing DSIPFT over bounding surface %s...",DSIMesh);
  else
   Log("Computing DSIPFT: (R,NPts)=(%e,%i)",DSIRadius, DSIPoints);

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(DSIMesh, DSIRadius, DSIPoints, false, GT);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT) GT->Apply(XTorque);

  double EpsAbs = TENTHIRDS / ZVAC;
  double  MuAbs = TENTHIRDS * ZVAC;

  /***************************************************************/
  /* get incident and scattered fields at the cubature points    */
  /***************************************************************/
  Log(" Computing incident fields at cubature points...");
  HMatrix *FInc  = G->GetFields(IF, 0, Omega, SCRMatrix);
  Log(" Computing scattered fields at cubature points...");
  HMatrix *FScat = G->GetFields( 0, J, Omega, SCRMatrix);

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  Log(" Evaluating cubature rule...");
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nr=0; nr<SCRMatrix->NR; nr++)
   { 
     double w, X[3], nHat[3];
     SCRMatrix->GetEntriesD(nr, "0:2", X);
     SCRMatrix->GetEntriesD(nr, "3:5", nHat);
     w = SCRMatrix->GetEntryD(nr, 6);

     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);

     cdouble ES[3], HS[3], ET[3], HT[3];
     FScat->GetEntries(nr, "0:2", ES);
     FScat->GetEntries(nr, "3:5", HS);
     FInc ->GetEntries(nr, "0:2", ET);
     FInc ->GetEntries(nr, "3:5", HT);
     for(int Mu=0; Mu<3; Mu++)
      { ET[Mu] += ES[Mu];
        HT[Mu] += HS[Mu];
      };

     // absorbed power 
     PFT[PFT_PABS]  -= 0.25 * w * (  HVMVP(ET, NMatrix[PFT_PABS], HT)
                                    -HVMVP(HT, NMatrix[PFT_PABS], ET)
                                  );
     // scattered power 
     PFT[PFT_PSCAT] += 0.25 * w * (  HVMVP(ES, NMatrix[PFT_PABS], HS)
                                    -HVMVP(HS, NMatrix[PFT_PABS], ES)
                                  );
     // force and torque
     for(int nq=PFT_XFORCE; nq<=PFT_ZTORQUE; nq++)
      PFT[nq] += 0.25 * w * ( EpsAbs*HVMVP(ET, NMatrix[nq], ET)
                             + MuAbs*HVMVP(HT, NMatrix[nq], HT)
                            );
   };
  Log(" Done!");

  delete FInc;
  delete FScat;
  delete SCRMatrix;

}

/***************************************************************/
/* Get power, force, and torque by the displaced               */
/* surface-integral method.                                    */
/***************************************************************/
#if 0
void GetDSIPFTTrace(SWGGeometry *G, HMatrix *R,
                    cdouble Omega, double PFT[NUMPFT],
                    char *DSIMesh, double DSIRadius, int DSIPoints,
                    GTransformation *GT)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (DSIMesh)
   Log("Computing DSIPFT trace over bounding surface %s...",DSIMesh);
  else
   Log("Computing DSIPFT trace: (R,NPts)=(%e,%i)",DSIRadius, DSIPoints);

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(DSIMesh, DSIRadius, DSIPoints, false, GT);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT) GT->Apply(XTorque);

  double EpsAbs = TENTHIRDS / ZVAC;
  double  MuAbs = TENTHIRDS * ZVAC;

  /***************************************************************/
  /* precompute 1BF fields at cubature points                    */
  /***************************************************************/
  Log(" Precomputing incident fields at cubature points...");
  HMatrix *FInc  = G->GetFields(IF, 0, Omega, SCRMatrix);
  Log(" Computing scattered fields at cubature points...");
  HMatrix *FScat = G->GetFields( 0, J, Omega, SCRMatrix);

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  Log(" Evaluating trace rule...");
#ifdef USE_OPENMP
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nbfa=0; nbfa<NBF; nbfa++)
   for(int nbfb=0; nbfb<NBF; nbfb++)
    for(int nr=0; nr<SCRMatrix->NR; nr++)
     { 
       double w, X[3], nHat[3];
       SCRMatrix->GetEntriesD(nr, "0:2", X);
       SCRMatrix->GetEntriesD(nr, "3:5", nHat);
       w = SCRMatrix->GetEntryD(nr, 6);

       double NMatrix[NUMPFT][3][3];
       GetNMatrices(nHat, X, XTorque, NMatrix);

       cdouble *EAlpha, *HAlpha, *EBeta, *HBeta;
       ES = FMatrix->GetEntries("0:2", ES);
       FScat->GetEntries(nr, "3:5", HS);
       FInc ->GetEntries(nr, "0:2", ET);
       FInc ->GetEntries(nr, "3:5", HT);
     for(int Mu=0; Mu<3; Mu++)
      { ET[Mu] += ES[Mu];
        HT[Mu] += HS[Mu];
      };

     // absorbed power 
     PFT[PFT_PABS]  -= 0.25 * w * (  HVMVP(ET, NMatrix[PFT_PABS], HT)
                                    -HVMVP(HT, NMatrix[PFT_PABS], ET)
                                  );
     // scattered power 
     PFT[PFT_PSCAT] += 0.25 * w * (  HVMVP(ES, NMatrix[PFT_PABS], HS)
                                    -HVMVP(HS, NMatrix[PFT_PABS], ES)
                                  );
     // force and torque
     for(int nq=PFT_XFORCE; nq<=PFT_ZTORQUE; nq++)
      PFT[nq] += 0.25 * w * ( EpsAbs*HVMVP(ET, NMatrix[nq], ET)
                             + MuAbs*HVMVP(HT, NMatrix[nq], HT)
                            );
   };
  Log(" Done!");

  delete FInc;
  delete FScat;
  delete SCRMatrix;

}
#endif

} // namespace buff
