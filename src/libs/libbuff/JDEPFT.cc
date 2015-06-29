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
 * JDEPFT.cc     -- libbuff class methods for computing power, force,
 *               -- and torque in classical deterministic scattering
 *               -- problems using the "J \dot E" formalism
 *
 * homer reid    -- 1/2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fenv.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0,1)

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace buff {

SWGVolume *ResolveNBF(SWGGeometry *G, int nbf, int *pno, int *pnf);
cdouble GetJJ(HVector *JVector, HMatrix *Rytov, int nbfa, int nbfb);

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PFTIntegrandData
 {
   cdouble k;
   bool NeedDerivatives;
   IncField *IF;
   double XTorque[3];

 } PFTIntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFInc(double *x, double *b, double Divb,
                        void *UserData, double *I)
{
  (void) Divb; // unused

  PFTIntegrandData *PFTIData=(PFTIntegrandData *)UserData;

  cdouble k    = PFTIData->k;
  double kabs  = abs(k);
  IncField *IF = PFTIData->IF; 

  // get fields and derivatives at eval point by finite-differencing
  cdouble EH[6], EHP[6], EHM[6], dEH[3][6];
  IF->GetFields(x, EH);
  for(int Mu=0; Mu<3; Mu++)
   { 
     double xTweaked[3];
     xTweaked[0]=x[0];
     xTweaked[1]=x[1];
     xTweaked[2]=x[2];

     double Delta = (x[Mu]==0.0 ? 1.0e-4 : 1.0e-4*fabs(x[Mu]));
     if ( kabs > 1.0 )
      Delta = fmin(Delta, 1.0e-4/kabs);

     xTweaked[Mu] += Delta;
     IF->GetFields(xTweaked, EHP);
     xTweaked[Mu] -= 2.0*Delta;
     IF->GetFields(xTweaked, EHM);

     for(int Nu=0; Nu<6; Nu++)
      dEH[Mu][Nu] = (EHP[Nu]-EHM[Nu])/(2.0*Delta);
   };

  double XmXT[3];
  VecSub(x, PFTIData->XTorque, XmXT);

  cdouble *zI = (cdouble *)I;
  memset(zI, 0, NUMPFTIS*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { zI[0] += b[Mu]*EH[Mu];
     zI[1] += b[Mu]*dEH[0][Mu];
     zI[2] += b[Mu]*dEH[1][Mu];
     zI[3] += b[Mu]*dEH[2][Mu];
     zI[4] += b[Mu]*(XmXT[1]*dEH[2][Mu]-XmXT[2]*dEH[1][Mu]);
     zI[5] += b[Mu]*(XmXT[2]*dEH[0][Mu]-XmXT[0]*dEH[2][Mu]);
     zI[6] += b[Mu]*(XmXT[0]*dEH[1][Mu]-XmXT[1]*dEH[0][Mu]);
   };

}

/***************************************************************/
/* compute PFT integrals between an SWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(SWGVolume *O, int nbf, IncField *IF,
                           cdouble Omega, cdouble IPFT[NUMPFTIS])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->IF = IF;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0;
  if (O->OTGT) O->OTGT->Apply(PFTIData->XTorque);
  if (O->GT) O->GT->Apply(PFTIData->XTorque);

  cdouble Error[NUMPFTIS];
  BFInt(O, nbf, PFTIntegrand_BFInc, (void *)PFTIData,
        2*NUMPFTIS, (double *)IPFT, (double *)Error,
        33, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetJDEPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HVector *RHSVector,
                   HMatrix *Rytov, HMatrix *PFTMatrix, bool *NeedFT)
{ 
  bool DefaultNeedFT[6]={true, true, true, true, true, true};
  if (NeedFT==0) NeedFT=DefaultNeedFT;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NO              = G->NumObjects;
  SWGVolume **Objects = G->Objects;
  int *BFIndexOffset  = G->BFIndexOffset;
  if (    (PFTMatrix==0)
       || (PFTMatrix->NR != NO)
       || (PFTMatrix->NC != NUMPFT)
     )
   ErrExit("invalid PFTMatrix in GetJDEPFT");

  /***************************************************************/
  /* precompute total power (extinction) if RHSVector is present */
  /***************************************************************/
  if (RHSVector)
   { cdouble PreFactor = -II*Omega*ZVAC;
     for(int no=0, nbf=0; no<NO; no++)
      for(int nf=0; nf<Objects[no]->NumInteriorFaces; nf++, nbf++)
       {
         cdouble EAlpha=PreFactor*RHSVector->GetEntry(nbf);
         cdouble jAlpha=JVector->GetEntry(nbf);
         PFTMatrix->AddEntry(no, 1, 0.5*real( conj(jAlpha) * EAlpha ) );
       };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumThreads=1;
#ifdef USE_OPENMP 
  NumThreads = GetNumThreads();
#endif
  int NQ = NUMPFT;
  int NONQ = NO*NQ;
  double *DeltaPFT = (double *)mallocEC(NumThreads*NONQ*sizeof(double));

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
  int TotalBFs = G->TotalBFs;
  int NumPairs = TotalBFs*(TotalBFs+1)/2;
  cdouble IKZ=II*Omega*ZVAC;
#ifdef USE_OPENMP
  Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      \
                         collapse(2),              \
                         num_threads(NumThreads)
#endif
  for(int nbfa=0; nbfa<TotalBFs; nbfa++) 
   for(int nbfb=0; nbfb<TotalBFs; nbfb++) 
    { 
      if (nbfb<nbfa) continue;
      if (nbfb==nbfa) LogPercent(nbfa*(nbfa+1)/2,NumPairs,10);

      int noa, nfa;
      SWGVolume *OA = ResolveNBF(G, nbfa, &noa, &nfa); 

      int nob, nfb;
      SWGVolume *OB = ResolveNBF(G, nbfb, &nob, &nfb);
   
      cdouble JJ = GetJJ(JVector, Rytov, nbfa, nbfb);
      if (JJ==0.0) continue;
   
      cdouble dG[6];
      FIBBICache *GCache=0, *dGCache=0;
      if (noa==nob)
       { GCache=G->ObjectGCaches[noa];
         dGCache=G->ObjectdGCaches[noa];
       };
      cdouble GG=GetGMatrixElement(OA, nfa, OB, nfb, Omega,
                                   GCache, dG, dGCache);

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      
       int Offset = nt*NONQ + noa*NQ;
       DeltaPFT[ Offset + 0 ] += 0.5*real( IKZ*JJ*GG );
       for(int Mu=0; Mu<6; Mu++)
        DeltaPFT[ Offset + 2 + Mu ] += 0.5*TENTHIRDS*imag( IKZ*JJ*dG[Mu] ) / real(Omega);

       if (nbfb>nbfa)
        { Offset = nt*NONQ + nob*NQ;
          DeltaPFT[ Offset + 0 ] += 0.5*real( IKZ*conj(JJ*GG) );
          for(int Mu=0; Mu<6; Mu++)
           DeltaPFT[ Offset + 2 + Mu ] += 0.5*TENTHIRDS*imag( IKZ*conj(JJ*dG[Mu]) ) / real(Omega);
        };

    }; // end of multithreaded loop
  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  for(int no=0; no<NO; no++)
   for(int nq=0; nq<NQ; nq++)
    for(int nt=0; nt<NumThreads; nt++)
     PFTMatrix->AddEntry(no, nq, DeltaPFT[ nt*NONQ + no*NQ + nq ]);

  /***************************************************************/
  /* add incident-field contributions ****************************/
  /***************************************************************/
  double Elapsed=Secs();
  if (IF)
   { 
      for(int no=0; no<NO; no++)
       { 
         SWGVolume *O = Objects[no];
         int Offset   = BFIndexOffset[no];
         int NBF      = O->NumInteriorFaces;

         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         double P=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Tx=0.0, Ty=0.0, Tz=0.0; 
#ifdef USE_OPENMP
NumThreads = GetNumThreads();
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:P, Fx, Fy, Fz, Tx, Ty, Tz)
#endif
         for(int nbf=0; nbf<NBF; nbf++)
          { 
            cdouble IPFT[7];
            GetPFTIntegrals_BFInc(O, nbf, IF, Omega, IPFT);
            cdouble J = conj(JVector->GetEntry(Offset + nbf));
            P  += real( J*IPFT[0] );
            Fx += imag( J*IPFT[1] );
            Fy += imag( J*IPFT[2] );
            Fz += imag( J*IPFT[3] );
            Tx += imag( J*IPFT[4] );
            Ty += imag( J*IPFT[5] );
            Tz += imag( J*IPFT[6] );
          };

         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         double Factor = 0.5;
         PFTMatrix->AddEntry(no, 0, Factor * P );
         Factor*=TENTHIRDS/real(Omega);
         PFTMatrix->AddEntry(no, 2, Factor * Fx );
         PFTMatrix->AddEntry(no, 3, Factor * Fy );
         PFTMatrix->AddEntry(no, 4, Factor * Fz );
         PFTMatrix->AddEntry(no, 5, Factor * Tx );
         PFTMatrix->AddEntry(no, 6, Factor * Ty );
         PFTMatrix->AddEntry(no, 7, Factor * Tz );
       };

   }; // if (IF)
  Elapsed = Secs() - Elapsed;
  //AddTaskTiming(5,Elapsed);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int no=0; no<NO; no++)
   PFTMatrix->AddEntry(no, 1, -1.0*PFTMatrix->GetEntry(no,0));

  return PFTMatrix;

}
  
} // namespace buff
