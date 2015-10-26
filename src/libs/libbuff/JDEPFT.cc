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
cdouble GetJJ(HVector *JVector, HMatrix *DMatrix, int nbfa, int nbfb);

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
  double k       = real(PFTIData->k);

  double XmXT[3]; 
  VecSub(xA, PFTIData->XTorque, XmXT);

  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB;
  double PEFIE         = DotProduct - ScalarProduct/(k*k);

  double R[3];
  VecSub(xA, xB, R);
  double r=VecNorm(R);
  if (r==0.0) 
   { 
     memset(I, 0, NUMPFT*sizeof(double));
     I[PFT_PABS] = PEFIE * k/(4.0*M_PI);
     return;
   };

  double kr  = k*r, coskr=cos(kr), sinkr=sin(kr);
  double ImG = sinkr/(4.0*M_PI*r);
  double f1  = (kr*coskr - sinkr) / (4.0*M_PI*r*r*r);

  I[PFT_PABS]    = PEFIE * ImG;
  I[PFT_XFORCE]  = PEFIE * R[0] * f1;
  I[PFT_YFORCE]  = PEFIE * R[1] * f1;
  I[PFT_ZFORCE]  = PEFIE * R[2] * f1;
  I[PFT_XTORQUE] = PEFIE * (XmXT[1]*R[2]-XmXT[2]*R[1]) * f1;
  I[PFT_YTORQUE] = PEFIE * (XmXT[2]*R[0]-XmXT[0]*R[2]) * f1;
  I[PFT_ZTORQUE] = PEFIE * (XmXT[0]*R[1]-XmXT[1]*R[0]) * f1;

  double k2  = k*k;
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
 
  I[PFT_XTORQUE] -=  DivbB*(bA[1]*R[2]-bA[2]*R[1])*f1/k2;
  I[PFT_YTORQUE] -=  DivbB*(bA[2]*R[0]-bA[0]*R[2])*f1/k2;
  I[PFT_ZTORQUE] -=  DivbB*(bA[0]*R[1]-bA[1]*R[0])*f1/k2;

  I[PFT_XTORQUE] += f2*bAxbB[0] + f3*bBdotR*bAxR[0];
  I[PFT_YTORQUE] += f2*bAxbB[1] + f3*bBdotR*bAxR[1];
  I[PFT_ZTORQUE] += f2*bAxbB[2] + f3*bBdotR*bAxR[2];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFInc(double *x, double *b, double Divb,
                        void *UserData, double I[NUMPFT])
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
  memset(zI, 0, NUMPFT*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { zI[PFT_PABS]    += b[Mu]*EH[Mu];
     zI[PFT_XFORCE]  += b[Mu]*dEH[0][Mu];
     zI[PFT_YFORCE]  += b[Mu]*dEH[1][Mu];
     zI[PFT_ZFORCE]  += b[Mu]*dEH[2][Mu];
     zI[PFT_XTORQUE] += b[Mu]*(XmXT[1]*dEH[2][Mu]-XmXT[2]*dEH[1][Mu]);
     zI[PFT_YTORQUE] += b[Mu]*(XmXT[2]*dEH[0][Mu]-XmXT[0]*dEH[2][Mu]);
     zI[PFT_ZTORQUE] += b[Mu]*(XmXT[0]*dEH[1][Mu]-XmXT[1]*dEH[0][Mu]);
   };

  zI[PFT_XTORQUE] += b[1]*EH[2] - b[2]*EH[1];
  zI[PFT_YTORQUE] += b[2]*EH[0] - b[0]*EH[2];
  zI[PFT_ZTORQUE] += b[0]*EH[1] - b[1]*EH[0];

}

/***************************************************************/
/* compute PFT integrals between an SWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(SWGVolume *O, int nbf, IncField *IF,
                           cdouble Omega, cdouble IBFInc[NUMPFT])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->IF = IF;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0;
  if (O->OTGT) O->OTGT->Apply(PFTIData->XTorque);
  if (O->GT) O->GT->Apply(PFTIData->XTorque);

  cdouble Error[NUMPFT];
  int fdim=2*NUMPFT;
  BFInt(O, nbf, PFTIntegrand_BFInc, (void *)PFTIData,
        fdim, (double *)IBFInc, (double *)Error, 33, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals_BFBF(SWGVolume *Oa, int nbfa,
                          SWGVolume *Ob, int nbfb,
                          cdouble Omega, double IBFBF[NUMPFT])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;

  double Error[NUMPFT];
  int ncv = CompareBFs(Oa, nbfa, Ob, nbfb);
  int NumPts = (ncv > 0) ? 33 : 16;
  int fdim = NUMPFT;
  BFBFInt(Oa, nbfa, Ob, nbfb, PFTIntegrand_BFBF, (void *)PFTIData,
          fdim, IBFBF, Error, NumPts, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddIFContributionsToJDEPFT(SWGGeometry *G, HVector *JVector,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix)
{
  if ( PFTMatrix->NR!=G->NumObjects || PFTMatrix->NC != NUMPFT )
   ErrExit("%s:%i: internal error", __FILE__, __LINE__);

  int NT=1;
#ifdef USE_OPENMP
  NT = GetNumThreads();
#endif
  int NO=G->NumObjects;
  int NQ=NUMPFT;
  double *PartialPFT=new double[NT*NO*NQ];
  memset(PartialPFT, 0, NT*NO*NQ*sizeof(double));

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int nbf=0; nbf<G->TotalBFs; nbf++)
   { 
     cdouble JStar = conj(JVector->GetEntry(nbf));

     int no, nf;
     SWGVolume *O = ResolveNBF(G, nbf, &no, &nf);

     cdouble IBFInc[NUMPFT];
     GetPFTIntegrals_BFInc(O, nf, IF, Omega, IBFInc);

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFT = PartialPFT + (nt*NO + no)*NUMPFT;
     dPFT[PFT_PABS] += real(JStar*IBFInc[PFT_PABS]);
     for(int Mu=0; Mu<6; Mu++)
      dPFT[PFT_XFORCE+Mu] += imag(JStar*IBFInc[PFT_XFORCE+Mu]);
   };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   double PFactor = 0.5;
   double FTFactor = 0.5*TENTHIRDS/real(Omega);
   for(int nt=0; nt<NT; nt++)
    for(int no=0; no<NO; no++)
     { 
       double *dPFT = PartialPFT + (nt*NO + no)*NUMPFT;
       PFTMatrix->AddEntry(no, PFT_PABS, PFactor*dPFT[PFT_PABS]);
       for(int nq=2; nq<NUMPFT; nq++)
        PFTMatrix->AddEntry(no, nq, FTFactor*dPFT[nq]);
     };

  delete[] PartialPFT;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetJDEPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HVector *RHSVector,
                   HMatrix *DMatrix, HMatrix *PFTMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NO = G->NumObjects;
  SWGVolume **Objects = G->Objects;
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
  if ( DeltaPFTSize < (NumThreads*NONQ) )
   { Log("(re)allocating DeltaPFT (%i,%i)",DeltaPFTSize,NumThreads*NONQ);
     DeltaPFTSize=NumThreads*NONQ;
     if (DeltaPFT) free(DeltaPFT);
     DeltaPFT = (double *)mallocEC(DeltaPFTSize*sizeof(double));
   };
  memset(DeltaPFT, 0, DeltaPFTSize*sizeof(double));

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
  int TotalBFs    = G->TotalBFs;
  double PPreFac  = real(Omega)*ZVAC;
  double FTPreFac = TENTHIRDS*ZVAC;
#ifdef USE_OPENMP
  Log("JDE OpenMP multithreading (%i threads)",NumThreads);
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
   
      cdouble JJ = GetJJ(JVector, DMatrix, nbfa, nbfb);
      if (JJ==0.0) continue;

      double IBFBF[NUMPFT];
      GetPFTIntegrals_BFBF(OA, nfa, OB, nfb, Omega, IBFBF);
      double ImG=IBFBF[0];
      double *ImdG=IBFBF+PFT_XFORCE;

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int OffsetA = nt*NONQ + noa*NQ;
      int OffsetB = nt*NONQ + nob*NQ;

       if (nbfa==nbfb)
        DeltaPFT[ OffsetA + PFT_PABS ] -= 0.5*PPreFac*real(JJ)*ImG;
       else if (noa==nob) // nbfb > nbfa but both on same object
        { 
          DeltaPFT[ OffsetA + PFT_PABS] -= PPreFac*real(JJ)*ImG;
          for(int Mu=0; Mu<6; Mu++)
           DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] -= FTPreFac*imag(JJ)*ImdG[Mu];
        }
       else // nbfb > nbfa and on different objects
        { 
          DeltaPFT[ OffsetA + PFT_PABS] -= 0.5*PPreFac*real(JJ)*ImG;
          DeltaPFT[ OffsetB + PFT_PABS] -= 0.5*PPreFac*real(JJ)*ImG;
          for(int Mu=0; Mu<6; Mu++)
           { DeltaPFT[ OffsetA + PFT_XFORCE + Mu ] -= 0.5*FTPreFac*imag(JJ)*ImdG[Mu];
             DeltaPFT[ OffsetB + PFT_XFORCE + Mu ] -= 0.5*FTPreFac*imag(JJ)*ImdG[Mu];
           };
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
  if (IF)
   AddIFContributionsToJDEPFT(G, JVector, IF, Omega, PFTMatrix);

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
