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

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace scuff {

HMatrix *GetSCRMatrix(char *BSMesh, double R, int NumPoints,
                      GTransformation *GT1, GTransformation *GT2);

void GetNMatrices(double nHat[3], double X[3], double XTorque[3],
                  double NMatrix[NUMPFT][3][3], 
                  bool *NeedQuantity=0);

double HVMVP(cdouble V1[3], double M[3][3], cdouble V2[3]);

                }

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace buff{

void Get1BFFields(SWGVolume *O, int nf, cdouble Omega, double X[3],
                  cdouble EH[6]);

cdouble HVMVP2(cdouble V1[3], double M[3][3], cdouble V2[3])
{ 
  cdouble Sum=0.0;

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    Sum += conj(V1[Mu]) * M[Mu][Nu] * V2[Nu];

  return Sum;
}

SWGVolume *ResolveNBF(SWGGeometry *G, int nbf, int *pno, int *pnf);

/***************************************************************/
/* Get power, force, and torque by the displaced               */
/* surface-integral method.                                    */
/***************************************************************/
void GetDSIPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
               HVector *JVector, double PFT[NUMPFT],
               GTransformation *GT1, GTransformation *GT2, 
               PFTOptions *Options)
{
  char *DSIMesh    = Options ? Options->DSIMesh    : 0;
  double DSIRadius = Options ? Options->DSIRadius  : 5.0;
  int DSIPoints    = Options ? Options->DSIPoints  : 110;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (DSIMesh)
   Log("DSI Computing DSIPFT over bounding surface %s...",DSIMesh);
  else
   Log("DSI Computing DSIPFT: (R,NPts)=(%e,%i)",DSIRadius, DSIPoints);

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(DSIMesh, DSIRadius, DSIPoints, GT1, GT2);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT1) GT1->Apply(XTorque);
  if (GT2) GT2->Apply(XTorque);

  double EpsAbs = TENTHIRDS / ZVAC;
  double  MuAbs = TENTHIRDS * ZVAC;

  /***************************************************************/
  /* get incident and scattered fields at the cubature points    */
  /***************************************************************/
  Log("DSI Computing incident fields at cubature points...");
  HMatrix *FInc  = IF ? G->GetFields(IF, 0, Omega, SCRMatrix) : 0;
  Log("DSI Computing scattered fields at cubature points...");
  HMatrix *FScat = JVector ? G->GetFields( 0, JVector, Omega, SCRMatrix) : 0;

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  Log("DSI Evaluating cubature rule...");
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nr=0; nr<SCRMatrix->NR; nr++)
   { 
     double w, X[3], nHat[3];
     SCRMatrix->GetEntriesD(nr, "0:2", X);
     SCRMatrix->GetEntriesD(nr, "3:5", nHat);
     w = SCRMatrix->GetEntryD(nr, 6);

     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);

     cdouble ES[3]={0.0,0.0,0.0};
     cdouble HS[3]={0.0,0.0,0.0};
     cdouble ET[3]={0.0,0.0,0.0};
     cdouble HT[3]={0.0,0.0,0.0};
     if (FScat)
      { FScat->GetEntries(nr, "0:2", ES);
        FScat->GetEntries(nr, "3:5", HS);
      };
     if (FInc)
      { FInc->GetEntries(nr, "0:2", ET);
        FInc->GetEntries(nr, "3:5", HT);
      };
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
  Log("DSI Done!");

  if (FInc) delete FInc;
  if (FScat) delete FScat;
  delete SCRMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetDSIPFTTrace(SWGGeometry *G, cdouble Omega, HMatrix *DRMatrix,
                    double PFT[NUMPFT],
                    GTransformation *GT1, GTransformation *GT2,
                    PFTOptions *Options)
{
  char *DSIMesh    = Options ? Options->DSIMesh    : 0;
  double DSIRadius = Options ? Options->DSIRadius  : 5.0;
  int DSIPoints    = Options ? Options->DSIPoints  : 110;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (DSIMesh)
   Log("DSI Computing DSIPFT trace over bounding surface %s...",DSIMesh);
  else
   Log("DSI Computing DSIPFT trace: (R,NPts)=(%e,%i)",DSIRadius, DSIPoints);

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(DSIMesh, DSIRadius, DSIPoints, GT1, GT2);

  /***************************************************************/
  /* precompute 1BF fields at cubature points                    */
  /***************************************************************/
  int NX  = SCRMatrix->NR;
  int NBF = G->TotalBFs;
  HMatrix *ehMatrix=new HMatrix(6, NX*NBF, LHM_COMPLEX);
  Log("DSI Precomputing 1BF fields at cubature points...");
  int NumThreads=1;
#ifdef USE_OPENMP
  NumThreads=GetNumThreads();
  Log("DSI using %i OpenMP threads",NumThreads);
#pragma omp parallel for num_threads(NumThreads), \
                         collapse(2),             \
                         schedule(dynamic,1)
#endif
  for(int nx=0; nx<NX; nx++)
   for(int nbf=0; nbf<NBF; nbf++)
    { 
      if (nbf==0 && G->LogLevel>=BUFF_VERBOSE_LOGGING) 
       LogPercent(nx,NX,10);

      double X[3];
      X[0] = SCRMatrix->GetEntryD(nx, 0);
      X[1] = SCRMatrix->GetEntryD(nx, 1);
      X[2] = SCRMatrix->GetEntryD(nx, 2);

      int no, nf;
      ResolveNBF(G, nbf, &no, &nf);
      cdouble *eh = ehMatrix->ZM + 6*(nx*NBF + nbf);
      Get1BFFields(G->Objects[no], nf, Omega, X, eh);
    };

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  Log("DSI Evaluating cubature rule...");
  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT1) GT1->Apply(XTorque);
  if (GT2) GT2->Apply(XTorque);
  double EpsAbs = TENTHIRDS / ZVAC;
  double  MuAbs = TENTHIRDS * ZVAC;
  cdouble *DeltaPFT=(cdouble *)mallocEC(NumThreads*NUMPFT*sizeof(cdouble));
#ifdef USE_OPENMP
  Log("DSI using %i OpenMP threads",NumThreads);
#pragma omp parallel for num_threads(NumThreads), \
                         schedule(dynamic,1)
#endif
  for(int nx=0; nx<NX; nx++)
   for(int nbfa=0; nbfa<NBF; nbfa++)
    { 
      if (nbfa==0) LogPercent(nx,NX,100);

      double X[3], nHat[3], w;
      X[0]    = SCRMatrix->GetEntryD(nx, 0);
      X[1]    = SCRMatrix->GetEntryD(nx, 1);
      X[2]    = SCRMatrix->GetEntryD(nx, 2);
      nHat[0] = SCRMatrix->GetEntryD(nx, 3);
      nHat[1] = SCRMatrix->GetEntryD(nx, 4);
      nHat[2] = SCRMatrix->GetEntryD(nx, 5);
      w       = SCRMatrix->GetEntryD(nx, 6);

      double NMatrix[NUMPFT][3][3];
      GetNMatrices(nHat, X, XTorque, NMatrix);

      cdouble *EA = ehMatrix->ZM + 6*(nx*NBF + nbfa);
      cdouble *HA = EA+3;

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      for(int nbfb=nbfa; nbfb<NBF; nbfb++)
       {
         cdouble *EB = ehMatrix->ZM + 6*(nx*NBF + nbfb);
         cdouble *HB = EB+3;

         cdouble JJ = DRMatrix->GetEntry(nbfb, nbfa);

         // absorbed/radiated power 
         DeltaPFT[nt*NUMPFT + PFT_PABS]
          -= 0.25 * w * JJ * (  HVMVP2(EA, NMatrix[PFT_PABS], HB)
                               -HVMVP2(HA, NMatrix[PFT_PABS], EB)
                             );

         // force/torque
         for(int nq=PFT_XFORCE; nq<=PFT_ZTORQUE; nq++)
          DeltaPFT[nt*NUMPFT + nq]
           += 0.25 * w * JJ * (  EpsAbs*HVMVP2(EA, NMatrix[nq], EB)
                                + MuAbs*HVMVP2(HA, NMatrix[nq], HB)
                              );

if (nbfb>nbfa)
{ JJ = conj(JJ);
         DeltaPFT[nt*NUMPFT + PFT_PABS]
          -= 0.25 * w * JJ * (  HVMVP2(EB, NMatrix[PFT_PABS], HA)
                               -HVMVP2(HB, NMatrix[PFT_PABS], EA)
                             );

         // force/torque
         for(int nq=PFT_XFORCE; nq<=PFT_ZTORQUE; nq++)
          DeltaPFT[nt*NUMPFT + nq]
           += 0.25 * w * JJ * (  EpsAbs*HVMVP2(EB, NMatrix[nq], EA)
                                + MuAbs*HVMVP2(HB, NMatrix[nq], HA)
                              );
};

       }; // for (nbfb=...)

    }; // for(int nx=0; ...), for(int nbfa=0; ...)

  Log("DSI Done!");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nq=0; nq<NUMPFT; nq++)
   for(int nt=0; nt<NumThreads; nt++)
    PFT[nq] += real(DeltaPFT[nt*NUMPFT + nq]);
   
  free(DeltaPFT);
  delete ehMatrix;
  delete SCRMatrix;

}

} // namespace buff{
