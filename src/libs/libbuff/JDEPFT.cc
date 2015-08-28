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
   IncField *IF;
   double XTorque[3];

 } PFTIntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFBF(double *xA, double *bA, double DivbA,
                       double *xB, double *bB, double DivbB,
                       void *UserData, double *I)
{
  (void) DivbA; // unused
  (void) DivbB; // unused

  PFTIntegrandData *PFTIData=(PFTIntegrandData *)UserData;

  double k    = real(PFTIData->k);

  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB;
  double PEFIE         = DotProduct - ScalarProduct/(k*k);

  double R[3];
  VecSub(xA, xB, R);
  double r=VecNorm(R);
  if (r==0.0) 
   { memset(I, 0, 6*sizeof(double));
     return;
   };

  double kr = k*r, coskr=cos(kr), sinkr=sin(kr);
  double f1 = (kr*coskr - sinkr) / (4.0*M_PI*r*r*r);

  I[0] = PEFIE * R[0] * f1;
  I[1] = PEFIE * R[1] * f1;
  I[2] = PEFIE * R[2] * f1;

  double kr2 = kr*kr;
  double f2  = (kr*coskr + (kr2-1.0)*sinkr) / (4.0*M_PI*kr2*r);
  double f3  = (-3.0*kr*coskr + (3.0-kr2)*sinkr) / (4.0*M_PI*kr2*r*r*r);

  double bAxbB[3], bAxR[3], bBdotR=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;
     bAxbB[Mu] = bA[MP1]*bB[MP2] - bA[MP2]*bB[MP1];
     bAxR[Mu]  = bA[MP1]*R[MP2]  - bA[MP2]*R[MP1];
     bBdotR   += bB[Mu]*R[Mu];
   };

  I[3] = f2*bAxbB[0] + f3*bBdotR*bAxR[0];
  I[4] = f2*bAxbB[1] + f3*bBdotR*bAxR[1];
  I[5] = f2*bAxbB[2] + f3*bBdotR*bAxR[2];

}

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
  memset(zI, 0, 7*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { zI[0] += b[Mu]*EH[Mu];
     zI[1] += b[Mu]*dEH[0][Mu];
     zI[2] += b[Mu]*dEH[1][Mu];
     zI[3] += b[Mu]*dEH[2][Mu];
     zI[4] += b[Mu]*(XmXT[1]*dEH[2][Mu]-XmXT[2]*dEH[1][Mu]);
     zI[5] += b[Mu]*(XmXT[2]*dEH[0][Mu]-XmXT[0]*dEH[2][Mu]);
     zI[6] += b[Mu]*(XmXT[0]*dEH[1][Mu]-XmXT[1]*dEH[0][Mu]);
   };

  zI[4 + 0] += b[1]*EH[2] - b[2]*EH[1];
  zI[4 + 1] += b[2]*EH[0] - b[0]*EH[2];
  zI[4 + 2] += b[0]*EH[1] - b[1]*EH[0];

}

/***************************************************************/
/* compute PFT integrals between an SWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(SWGVolume *O, int nbf, IncField *IF,
                           cdouble Omega, cdouble IPFT[7])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->IF = IF;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0;
  if (O->OTGT) O->OTGT->Apply(PFTIData->XTorque);
  if (O->GT) O->GT->Apply(PFTIData->XTorque);

  cdouble Error[7];
  BFInt(O, nbf, PFTIntegrand_BFInc, (void *)PFTIData,
        14, (double *)IPFT, (double *)Error, 33, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals_BFBF(SWGVolume *Oa, int nbfa,
                          SWGVolume *Ob, int nbfb,
                          cdouble Omega, double IPFT[6])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;

  double Error[6];
  int ncv = CompareBFs(Oa, nbfa, Ob, nbfb);
  int NumPts = (ncv > 0 ) ? 33 : 16;
  BFBFInt(Oa, nbfa, Ob, nbfb, PFTIntegrand_BFBF, (void *)PFTIData,
          6, IPFT, Error, NumPts, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetJDEPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HVector *RHSVector,
                   HMatrix *Rytov, HMatrix *PFTMatrix)
{ 
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

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumThreads=1;
#ifdef USE_OPENMP 
  NumThreads = GetNumThreads();
#endif
  int NQ = NUMPFT;
  int NONQ = NO*NQ;
  static int DeltaPFTSize=0;
  static double *DeltaPFT=0;
  if ( DeltaPFTSize != (NumThreads*NONQ) )
   { DeltaPFTSize=NumThreads*NONQ;
     if (DeltaPFT) free(DeltaPFT);
     DeltaPFT = (double *)mallocEC(DeltaPFTSize*sizeof(double));
   };

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
  int TotalBFs    = G->TotalBFs;
  double PPreFac  = real(Omega)*ZVAC;
  double FTPreFac = TENTHIRDS*ZVAC;
#ifdef USE_OPENMP
  Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      	\
                         num_threads(NumThreads)
#endif
  for(int nbfa=0; nbfa<TotalBFs; nbfa++)
   for(int nbfb=nbfa; nbfb<TotalBFs; nbfb++)
    { 
      //if (nbfb==nbfa) LogPercent(nbfa*(nbfa+1)/2,NumPairs,100);
      if (nbfb==nbfa) LogPercent(nbfa, TotalBFs, 10);

      int noa, nfa;
      SWGVolume *OA = ResolveNBF(G, nbfa, &noa, &nfa);

      int nob, nfb;
      SWGVolume *OB = ResolveNBF(G, nbfb, &nob, &nfb);
   
      cdouble JJ = GetJJ(JVector, Rytov, nbfa, nbfb);
      if (JJ==0.0) continue;
   
      FIBBICache *GCache = (noa==nob) ? G->ObjectGCaches[noa] : 0;
      cdouble GG=GetGMatrixElement(OA, nfa, OB, nfb, Omega, GCache);

      double ImdG[6];
      GetPFTIntegrals_BFBF(OA, nfa, OB, nfb, Omega, ImdG);

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int Offset = nt*NONQ + noa*NQ;

       if (nbfa==nbfb)
        DeltaPFT[ Offset + PFT_PABS ] -= 0.5*PPreFac*real(JJ)*imag(GG);
       else // nbfb > nbfa
        { DeltaPFT[ Offset + PFT_PABS] -= PPreFac*real(JJ)*imag(GG);
          for(int Mu=0; Mu<6; Mu++)
           DeltaPFT[ Offset + PFT_XFORCE + Mu ] -= FTPreFac*imag(JJ)*ImdG[Mu];
        };

    }; // end of multithreaded loop
  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  PFTMatrix->Zero();
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
         PFTMatrix->AddEntry(no, PFT_PABS, Factor * P );
         Factor*=TENTHIRDS/real(Omega);
         PFTMatrix->AddEntry(no, PFT_XFORCE, Factor * Fx );
         PFTMatrix->AddEntry(no, PFT_YFORCE, Factor * Fy );
         PFTMatrix->AddEntry(no, PFT_ZFORCE, Factor * Fz );
         PFTMatrix->AddEntry(no, PFT_XTORQUE, Factor * Tx );
         PFTMatrix->AddEntry(no, PFT_YTORQUE, Factor * Ty );
         PFTMatrix->AddEntry(no, PFT_ZTORQUE, Factor * Tz );
       };

   }; // if (IF)
  Elapsed = Secs() - Elapsed;
  //AddTaskTiming(5,Elapsed);

  /***************************************************************/
  /* compute scattered power only if RHSVector is present        */
  /* PScat = PTotal - PAbsorbed                                  */
  /***************************************************************/
  if (RHSVector)
   { 
     cdouble PreFactor = -II*Omega*ZVAC;
     for(int no=0, nbf=0; no<NO; no++)
      { 
        PFTMatrix->SetEntry(no, PFT_PSCAT, -1.0*PFTMatrix->GetEntry(no,PFT_PABS));
        for(int nf=0; nf<Objects[no]->NumInteriorFaces; nf++, nbf++)
         { cdouble EAlpha=PreFactor*RHSVector->GetEntry(nbf);
           cdouble jAlpha=JVector->GetEntry(nbf);
           PFTMatrix->AddEntry(no, PFT_PSCAT, 0.5*real( conj(jAlpha) * EAlpha ) );
         };
      };
   };

  return PFTMatrix;

}
  
} // namespace buff
