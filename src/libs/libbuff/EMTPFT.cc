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
 * EMTPFT.cc     -- libbuff class methods for computing power, force,
 *               -- and torque in classical deterministic scattering
 *               -- problems using the "energy/momentum transfer" approach
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

#define PFT_XTORQUE1 5
#define PFT_YTORQUE1 6
#define PFT_ZTORQUE1 7
#define PFT_XTORQUE2 8
#define PFT_YTORQUE2 9 
#define PFT_ZTORQUE2 10

#define NUMPFTT 11
#define NUMFTT  9

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
namespace buff {

SWGVolume *ResolveNBF(SWGGeometry *G, int nbf, int *pno, int *pnf);
cdouble GetJJ(HVector *JVector, HMatrix *DRMatrix, int nbfa, int nbfb);

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

  double XT[3]; 
  VecSub(xA, PFTIData->XTorque, XT);

  double R[3];
  VecSub(xA, xB, R);
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r=sqrt(r2), kr=k*r, kr2 = kr*kr, k2=k*k;

  /***************************************************************/
  /* polynomial factors ******************************************/
  /***************************************************************/
  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB;
  double PEFIE         = DotProduct - ScalarProduct/(k*k);
  double bAxbB[3], bAxR[3], XTxR[3], bBdR=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;
     bAxbB[Mu] = bA[MP1]*bB[MP2] - bA[MP2]*bB[MP1];
     bAxR[Mu]  = bA[MP1]*R[MP2]  - bA[MP2]*R[MP1];
     XTxR[Mu]  = XT[MP1]*R[MP2]  - XT[MP2]*R[MP1];
     bBdR     += bB[Mu]*R[Mu];
   };

  /***************************************************************/
  /* kernel factors **********************************************/
  /***************************************************************/
  double ImPhi, ImPsi, ImZeta;
  if ( fabs(kr) < 1.0e-3 )
   { double k3=k2*k, k5=k3*k2;
     ImPhi  =  (1.0 - kr2/6.0)  * k/(4.0*M_PI);
     ImPsi  = -(1.0 - kr2/10.0) * k3/(12.0*M_PI);
     ImZeta =  (1.0 - kr2/14.0) * k5/(60.0*M_PI);
   }
  else
   { double CosKR = cos(kr), SinKR=sin(kr), r3=r*r2, r5=r3*r2;
     ImPhi  = SinKR/(4.0*M_PI*r);
     ImPsi  = (kr*CosKR - SinKR)/(4.0*M_PI*r3);
     ImZeta = (-3.0*kr*CosKR + (3.0-kr2)*SinKR)/(4.0*M_PI*r5);
   };
  double PPPoK2 = ImPhi + ImPsi/k2;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  I[PFT_PABS] += PEFIE * ImPhi;
  for(int Mu=0; Mu<3; Mu++)
   { I[PFT_XFORCE   + Mu] += PEFIE * R[Mu] * ImPsi;
     I[PFT_XTORQUE1 + Mu] += bAxbB[Mu]*PPPoK2 + bAxR[Mu]*bBdR*ImZeta/k2;
     I[PFT_XTORQUE2 + Mu] += PEFIE * XTxR[Mu] * ImPsi;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFInc(double *x, double *b, double Divb,
                        void *UserData, double I[NUMPFTT])
{
  (void) Divb; // unused

  PFTIntegrandData *PFTIData=(PFTIntegrandData *)UserData;

  IncField *IF = PFTIData->IF;

  double XT[3];
  VecSub(x, PFTIData->XTorque, XT);

  // get fields and derivatives at eval point
  cdouble EH[6], dEH[3][6];
  IF->GetFields(x, EH);
  IF->GetFieldGradients(x, dEH);

  cdouble *zI = (cdouble *)I;
  memset(zI, 0, NUMPFT*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { zI[PFT_PABS]     += b[Mu]*EH[Mu];
     zI[PFT_XFORCE]   += b[Mu]*dEH[0][Mu];
     zI[PFT_YFORCE]   += b[Mu]*dEH[1][Mu];
     zI[PFT_ZFORCE]   += b[Mu]*dEH[2][Mu];
     zI[PFT_XTORQUE2] += b[Mu]*(XT[1]*dEH[2][Mu]-XT[2]*dEH[1][Mu]);
     zI[PFT_YTORQUE2] += b[Mu]*(XT[2]*dEH[0][Mu]-XT[0]*dEH[2][Mu]);
     zI[PFT_ZTORQUE2] += b[Mu]*(XT[0]*dEH[1][Mu]-XT[1]*dEH[0][Mu]);
   };

  zI[PFT_XTORQUE1] = b[1]*EH[2] - b[2]*EH[1];
  zI[PFT_YTORQUE1] = b[2]*EH[0] - b[0]*EH[2];
  zI[PFT_ZTORQUE1] = b[0]*EH[1] - b[1]*EH[0];

}

/***************************************************************/
/* compute PFT integrals between an SWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetPFTIntegrals_BFInc(SWGVolume *O, int nbf, IncField *IF,
                           cdouble Omega, cdouble IBFInc[NUMPFTT])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;
  PFTIData->IF = IF;
  PFTIData->XTorque[0]=PFTIData->XTorque[1]=PFTIData->XTorque[2]=0.0;
  if (O->OTGT) O->OTGT->Apply(PFTIData->XTorque);
  if (O->GT) O->GT->Apply(PFTIData->XTorque);

  cdouble Error[NUMPFTT];
  int IDim=2*NUMPFTT;
  int NumPts=33;
  BFInt(O, nbf, PFTIntegrand_BFInc, (void *)PFTIData,
        IDim, (double *)IBFInc, (double *)Error, NumPts, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPFTIntegrals_BFBF(SWGVolume *Oa, int nbfa,
                          SWGVolume *Ob, int nbfb,
                          cdouble Omega, double IBFBF[NUMPFTT])
{
  PFTIntegrandData MyPFTIData, *PFTIData=&MyPFTIData;
  PFTIData->k  = Omega;

  double Error[NUMPFTT];
  int ncv = CompareBFs(Oa, nbfa, Ob, nbfb);
  int NumPts = (ncv > 0) ? 33 : 16;
  int IDim = NUMPFTT;
  BFBFInt(Oa, nbfa, Ob, nbfb, PFTIntegrand_BFBF, (void *)PFTIData,
          IDim, IBFBF, Error, NumPts, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddIFContributionsToEMTPFTT(SWGGeometry *G, HVector *JVector,
                                 IncField *IF, cdouble Omega,
                                 HMatrix *PFTTMatrix)
{
  if ( PFTTMatrix->NR!=G->NumObjects || PFTTMatrix->NC != NUMPFT )
   ErrExit("%s:%i: internal error", __FILE__, __LINE__);

  int NT=1;
#ifdef USE_OPENMP
  NT = GetNumThreads();
#endif
  int NO=G->NumObjects;
  int NQ=NUMPFTT;
  double *PartialPFTT=new double[NT*NO*NQ];
  memset(PartialPFTT, 0, NT*NO*NQ*sizeof(double));

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int nbf=0; nbf<G->TotalBFs; nbf++)
   { 
     cdouble JStar = conj(JVector->GetEntry(nbf));

     int no, nf;
     SWGVolume *O = ResolveNBF(G, nbf, &no, &nf);

     cdouble IBFInc[NUMPFTT];
     GetPFTIntegrals_BFInc(O, nf, IF, Omega, IBFInc);

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFTT = PartialPFTT + (nt*NO + no)*NUMPFTT;
     dPFTT[PFT_PABS] += real(JStar*IBFInc[PFT_PABS]);
     for(int Mu=0; Mu<NUMFTT; Mu++)
      dPFTT[PFT_XFORCE+Mu] += imag(JStar*IBFInc[PFT_XFORCE+Mu]);
   };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   double PFactor = 0.5;
   double FTFactor = 0.5*TENTHIRDS/real(Omega);
   for(int nt=0; nt<NT; nt++)
    for(int no=0; no<NO; no++)
     { 
       double *dPFTT = PartialPFTT + (nt*NO + no)*NUMPFTT;
       PFTTMatrix->AddEntry(no, PFT_PABS, PFactor*dPFTT[PFT_PABS]);
       for(int nq=2; nq<NUMPFTT; nq++)
        PFTTMatrix->AddEntry(no, nq, FTFactor*dPFTT[nq]);
     };

  delete[] PartialPFTT;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddIFContributionsToEMTPFT(SWGGeometry *G, HVector *JVector,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix)
{
  int NO = PFTMatrix->NR;

  static HMatrix *PFTTMatrix=0;
  if (PFTTMatrix==0 || PFTTMatrix->NR != NO )
   { if (PFTTMatrix) delete PFTTMatrix;
     PFTTMatrix=new HMatrix(NO, NUMPFTT);
   };

  AddIFContributionsToEMTPFTT(G, JVector, IF, Omega, PFTTMatrix);
  for(int no=0; no<NO; no++)
   for(int nq=0; nq<NUMPFT; nq++)
    { PFTMatrix->SetEntry(no, nq, PFTTMatrix->GetEntry(no, nq));
      if (PFT_XTORQUE<=nq && nq<=PFT_ZTORQUE)
       PFTMatrix->AddEntry(no, nq, PFTTMatrix->GetEntry(no, nq+3));
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetEMTPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HVector *RHSVector,
                   HMatrix *DRMatrix, HMatrix *PFTMatrix)
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
   ErrExit("invalid PFTMatrix in GetEMTPFT");

  /***************************************************************/
  /* PFTTByObject[no] = contributions of object #no to PFTT      */
  /***************************************************************/
  static int NOSave=0;
  static HMatrix **PFTTByObject=0, *IncidentPFTT=0;
  if (NOSave!=NO)
   { if (PFTTByObject)
      { for(int no=0; no<NOSave; no++)
         if (PFTTByObject[no]) 
          delete PFTTByObject[no];
        free(PFTTByObject);
       if (IncidentPFTT)
        delete IncidentPFTT;
      };
     NOSave=NO;
     PFTTByObject=(HMatrix **)mallocEC(NO*sizeof(HMatrix));
     for(int no=0; no<NO; no++)
      PFTTByObject[no]=new HMatrix(NO, NUMPFTT);
     IncidentPFTT=new HMatrix(NO, NUMPFTT);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NT=1;
#ifdef USE_OPENMP 
  NT= GetNumThreads();
#endif
  int NQ      = NUMPFTT;
  int NONQ    = NO*NQ;
  int NO2NQ   = NO*NONQ;
  int NTNO2NQ = NT*NO2NQ;
  static int DeltaPFTTSize=0;
  static double *DeltaPFTT=0;
  if ( DeltaPFTTSize < NTNO2NQ )
   { Log("(re)allocating DeltaPFTT (%i,%i)",DeltaPFTTSize,NTNO2NQ);
     DeltaPFTTSize=NTNO2NQ;
     if (DeltaPFTT) free(DeltaPFTT);
     DeltaPFTT = (double *)mallocEC(DeltaPFTTSize*sizeof(double));
   };
  memset(DeltaPFTT, 0, DeltaPFTTSize*sizeof(double));

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
  int TotalBFs    = G->TotalBFs;
  double PPreFac  = real(Omega)*ZVAC;
  double FTPreFac = TENTHIRDS*ZVAC;
  bool UseSymmetry=true;
  char *s=getenv("BUFF_EMTPFT_SYMMETRY");
  if (s && s[0]=='0')
   { UseSymmetry=false;
     Log("Not using symmetry in EMTPFT calculation.");
   };
#ifdef USE_OPENMP
  Log("EMT OpenMP multithreading (%i threads)",NT);
#pragma omp parallel for schedule(dynamic,1),      	\
                         num_threads(NT)
#endif
  for(int nbfa=0; nbfa<TotalBFs; nbfa++)
   for(int nbfb=(UseSymmetry ? nbfa : 0); nbfb<TotalBFs; nbfb++)
    { 
      if (nbfb==(UseSymmetry ? nbfa : 0)) 
       LogPercent(nbfa, TotalBFs, 10);

      int noa, nfa;
      SWGVolume *OA = ResolveNBF(G, nbfa, &noa, &nfa);

      int nob, nfb;
      SWGVolume *OB = ResolveNBF(G, nbfb, &nob, &nfb);
   
      cdouble JJ = GetJJ(JVector, DRMatrix, nbfa, nbfb);
      if (JJ==0.0) continue;

      double IBFBF[NUMPFTT];
      GetPFTIntegrals_BFBF(OA, nfa, OB, nfb, Omega, IBFBF);
      double QPAbs=IBFBF[PFT_PABS];
      double *QFTT=IBFBF+PFT_XFORCE;

      double dPFTT[NUMPFTT];
      dPFTT[PFT_PABS] = -0.5*PPreFac*real(JJ)*QPAbs;
      for(int Mu=0; Mu<NUMFTT; Mu++)
       dPFTT[PFT_XFORCE + Mu ] -= 0.5*FTPreFac*imag(JJ)*QFTT[Mu];

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int OffsetA = nt*NO2NQ + noa*NONQ + nob*NQ;
      int OffsetB = nt*NO2NQ + nob*NONQ + noa*NQ; 

      if (!UseSymmetry)
       VecPlusEquals(DeltaPFTT + OffsetA, 1.0, dPFTT, NUMPFTT);
      else
       { if (nbfa==nbfb)
          DeltaPFTT[ OffsetA + PFT_PABS ] += dPFTT[PFT_PABS];
         else if (noa==nob) // nbfb > nbfa but both on same object
          VecPlusEquals(DeltaPFTT+OffsetA, 2.0, dPFTT, NUMPFTT);
         else // nbfb > nbfa and on different objects
          { 
            DeltaPFTT[OffsetA + PFT_PABS] += dPFTT[PFT_PABS];
            DeltaPFTT[OffsetB + PFT_PABS] += dPFTT[PFT_PABS];
            for(int Mu=0; Mu<NUMFTT; Mu++)
             { DeltaPFTT[ OffsetA + PFT_XFORCE + Mu ] += dPFTT[PFT_XFORCE + Mu];
               DeltaPFTT[ OffsetB + PFT_XFORCE + Mu ] -= dPFTT[PFT_XFORCE + Mu];
             };
          };
       };

    }; // end of multithreaded loop

  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  for(int no=0; no<NO; no++)
   PFTTByObject[no]->Zero();
  for(int noa=0; noa<NO; noa++)
   for(int nob=0; nob<NO; nob++)
    for(int nq=0; nq<NQ; nq++)
     for(int nt=0; nt<NT; nt++)
      PFTTByObject[nob]->AddEntry(noa, nq, DeltaPFTT[ nt*NO2NQ + noa*NONQ + nob*NQ + nq ]);

  /***************************************************************/
  /* get incident-field contributions ****************************/
  /***************************************************************/
  IncidentPFTT->Zero();
  if (IF)
   AddIFContributionsToEMTPFTT(G, JVector, IF, Omega, IncidentPFTT);
   
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTMatrix->Zero();
  for(int noa=0; noa<NO; noa++)
   for(int nq=0; nq<NUMPFT; nq++)
    { 
      for(int nob=0; nob<NO; nob++)
       PFTMatrix->AddEntry(noa,nq,PFTTByObject[nob]->GetEntry(noa,nq));

      PFTMatrix->AddEntry(noa,nq,IncidentPFTT->GetEntry(noa,nq));

      if (PFT_XTORQUE<=nq && nq<=PFT_ZTORQUE)
       { for(int nob=0; nob<NO; nob++)
          PFTMatrix->AddEntry(noa, nq, PFTTByObject[nob]->GetEntry(noa,nq+3));
         PFTMatrix->AddEntry(noa,nq,IncidentPFTT->GetEntry(noa,nq+3));
       };
    };

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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *ss=getenv("BUFF_ITEMIZE_PFT");
  if (ss && ss[0]=='1')
   { 
     static bool WrotePreamble=false;
     for(int noa=0; noa<NO; noa++)
      { FILE *f=vfopen("%s.%s.EMTPFT","a",
                        GetFileBase(G->GeoFileName),
                        G->Objects[noa]->Label);
        if (!f) continue;
        if (!WrotePreamble)
         { fprintf(f,"# EMTPFT contributions to object %s\n",G->Objects[noa]->Label);
           fprintf(f,"# columns: \n");
           fprintf(f,"# 1 frequency \n");
           fprintf(f,"# 2 destination object label \n");
           fprintf(f,"# 03-10 PAbs, PScat, Fxyz, Txyz (total)\n");
           fprintf(f,"# 11-21 PAbs, PScat, Fxyz, T1xyz, T2xyz (extinction)\n");
           int nc=22;
           for(int nob=0; nob<NO; nob++, nc+=NUMPFTT)
            fprintf(f,"# %i-%i PAbs, PScat, Fxyz, Txyz (object %s)\n",nc,nc+NUMPFTT-1,G->Objects[nob]->Label);
         };
        fprintf(f,"%e %s ",real(Omega),G->Objects[noa]->Label);
        for(int nq=0; nq<NUMPFT; nq++)
         fprintf(f,"%e ",PFTMatrix->GetEntryD(noa,nq));
        for(int nq=0; nq<NUMPFTT; nq++)
         fprintf(f,"%e ",IncidentPFTT->GetEntryD(noa,nq));
        for(int nob=0; nob<NO; nob++)
         for(int nq=0; nq<NUMPFTT; nq++)
          fprintf(f,"%e ",PFTTByObject[nob]->GetEntryD(noa,nq));
        fprintf(f,"\n");
        fclose(f);
      };
     WrotePreamble=true;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return PFTMatrix;

}
  
} // namespace buff
