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
typedef struct PFTIData
 {
   double Omega;
   IncField *IF;
   double *TorqueCenterA;
   double *TorqueCenterB;
   bool SameObject;

 } PFTIData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ExtinctionPFTIntegrand(double *x, double *b, double Divb,
                            void *UserData, double *I)
{
  (void) Divb; // unused

  PFTIData *Data       = (PFTIData *)UserData;
  IncField *IF         = Data->IF;
  double *TorqueCenter = Data->TorqueCenterA;

  double XT[3];
  VecSub(x, TorqueCenter, XT);

  // get fields and derivatives at eval point
  cdouble EH[6], dEH[3][6];
  IF->GetFields(x, EH);
  IF->GetFieldGradients(x, dEH);

  cdouble *Q = (cdouble *)I;
  memset(Q, 0, NUMPFTT*sizeof(cdouble));
  for(int Mu=0; Mu<3; Mu++)
   { Q[PFT_PABS]     += b[Mu]*EH[Mu];
     Q[PFT_XFORCE]   += b[Mu]*dEH[0][Mu];
     Q[PFT_YFORCE]   += b[Mu]*dEH[1][Mu];
     Q[PFT_ZFORCE]   += b[Mu]*dEH[2][Mu];
     Q[PFT_XTORQUE2] += b[Mu]*(XT[1]*dEH[2][Mu]-XT[2]*dEH[1][Mu]);
     Q[PFT_YTORQUE2] += b[Mu]*(XT[2]*dEH[0][Mu]-XT[0]*dEH[2][Mu]);
     Q[PFT_ZTORQUE2] += b[Mu]*(XT[0]*dEH[1][Mu]-XT[1]*dEH[0][Mu]);
   };

  Q[PFT_XTORQUE1] = b[1]*EH[2] - b[2]*EH[1];
  Q[PFT_YTORQUE1] = b[2]*EH[0] - b[0]*EH[2];
  Q[PFT_ZTORQUE1] = b[0]*EH[1] - b[1]*EH[0];

}

/***************************************************************/
/* compute PFT integrals between an SWG basis function and an  */
/* external field                                              */
/***************************************************************/
void GetExtinctionPFTIntegrals(SWGVolume *O, int nbf, IncField *IF,
                               cdouble Omega, cdouble Q[NUMPFTT])
{
  PFTIData MyData, *Data=&MyData;
  Data->Omega         = real(Omega);
  Data->IF            = IF;
  Data->TorqueCenterA = O->Origin;

  double Error[2*NUMPFTT];
  int IDim=2*NUMPFTT;
  int NumPts=33;
  BFInt(O, nbf, ExtinctionPFTIntegrand, (void *)Data,
        IDim, (double *)Q, Error, NumPts, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetExtinctionPFTT(SWGGeometry *G, HVector *JVector,
                       IncField *IF, cdouble Omega,
                       HMatrix *PFTTMatrix)
{
  if ( PFTTMatrix->NR!=G->NumObjects || PFTTMatrix->NC != NUMPFTT )
   ErrExit("%s:%i: internal error", __FILE__, __LINE__);

  int NT=1;
#ifdef USE_OPENMP
  NT = GetNumThreads();
#endif
  int NO=G->NumObjects;
  int NQ=NUMPFTT;
  int NTNONQ=NT*NO*NQ;
  
  static int DeltaPFTTSize=0;
  static double *DeltaPFTT=0;
  if (DeltaPFTTSize < NTNONQ)
   { DeltaPFTTSize=NTNONQ;
     DeltaPFTT=(double *)reallocEC(DeltaPFTT, DeltaPFTTSize*sizeof(double));
   };
  memset(DeltaPFTT, 0, NTNONQ*sizeof(double));

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
  for(int nbf=0; nbf<G->TotalBFs; nbf++)
   { 
     cdouble JStar = conj(JVector->GetEntry(nbf));

     int no, nf;
     SWGVolume *O = ResolveNBF(G, nbf, &no, &nf);

     cdouble Q[NUMPFTT];
     GetExtinctionPFTIntegrals(O, nf, IF, Omega, Q);

     int nt=0;
#ifdef USE_OPENMP
     nt = omp_get_thread_num();
#endif
     double *dPFTT = DeltaPFTT + (nt*NO + no)*NUMPFTT;
     dPFTT[PFT_PABS] += real(JStar*Q[PFT_PABS]);
     for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
      dPFTT[nq] += imag(JStar*Q[nq]);
   };

   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/
   double PFactor = 0.5;
   double FTFactor = 0.5*TENTHIRDS/real(Omega);
   PFTTMatrix->Zero();
   for(int nt=0; nt<NT; nt++)
    for(int no=0; no<NO; no++)
     { 
       double *dPFTT = DeltaPFTT + (nt*NO + no)*NUMPFTT;
       PFTTMatrix->AddEntry(no, PFT_PABS, PFactor*dPFTT[PFT_PABS]);
       for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
        PFTTMatrix->AddEntry(no, nq, FTFactor*dPFTT[nq]);
     };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ScatteredPFTIntegrand(double *xA, double *bA, double DivbA,
                           double *xB, double *bB, double DivbB,
                           void *UserData, double *I)
{
  (void) DivbA; // unused
  (void) DivbB; // unused

  PFTIData *Data        = (PFTIData*)UserData;
  double Omega          = Data->Omega;
  double k              = Data->Omega;
  bool SameObject       = Data->SameObject;
  double *TorqueCenterA = Data->TorqueCenterA;
  double *TorqueCenterB = Data->TorqueCenterB;

  double XTA[3], XTB[3];
  VecSub(xA, TorqueCenterA, XTA);
  VecSub(xB, TorqueCenterB, XTB);

  double R[3];
  VecSub(xA, xB, R);
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r=sqrt(r2), kr=k*r, kr2 = kr*kr, k2=k*k;

  /***************************************************************/
  /* polynomial factors ******************************************/
  /***************************************************************/
  double DotProduct    = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  double ScalarProduct = DivbA*DivbB;
  double PEFIE         = DotProduct - ScalarProduct/k2;
  double bAxbB[3], XTAxR[3], XTBxR[3];
  for(int Mu=0; Mu<3; Mu++)
   { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;
     bAxbB[Mu] = bA[MP1]*bB[MP2] - bA[MP2]*bB[MP1];
     XTAxR[Mu]  =  XTA[MP1]*R[MP2]  - XTA[MP2]*R[MP1];
     XTBxR[Mu]  = -XTB[MP1]*R[MP2]  + XTB[MP2]*R[MP1];
   };

  /***************************************************************/
  /* kernel factors **********************************************/
  /***************************************************************/
  cdouble Phi, Psi;
  if (SameObject)
   { double ImPhi, ImPsi;
     if ( fabs(kr) < 1.0e-3 )
      { double k3=k2*k;
        ImPhi  =  (1.0 - kr2/6.0)  * k/(4.0*M_PI);
        ImPsi  = -(1.0 - kr2/10.0) * k3/(12.0*M_PI);
      }
     else
      { double CosKR = cos(kr), SinKR=sin(kr), r3=r*r2;
        ImPhi  = SinKR/(4.0*M_PI*r);
        ImPsi  = (kr*CosKR - SinKR)/(4.0*M_PI*r3);
      };
     Phi=II*ImPhi;
     Psi=II*ImPsi;
   }
  else
   { cdouble ikr=II*kr;
     Phi  = exp(ikr)/(4.0*M_PI*r);
     Psi  = Phi*(ikr-1.0)/r2;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble *Q=(cdouble *)I;
  Q[PFT_PABS ] = 0.0;
  Q[PFT_PSCAT] = Omega * PEFIE * Phi;
  for(int Mu=0; Mu<3; Mu++)
   { Q[PFT_XFORCE   + Mu] = TENTHIRDS*PEFIE * R[Mu] * Psi;
     Q[PFT_XTORQUE1 + Mu] = TENTHIRDS*bAxbB[Mu]*Phi; //PPPoK2 + bAxR[Mu]*bBdR*ImZeta/k2;
     Q[PFT_XTORQUE2 + Mu] = TENTHIRDS*PEFIE * XTAxR[Mu] * Psi;
     Q[PFT_XTORQUE2 + 3 * Mu] = TENTHIRDS*PEFIE * XTBxR[Mu] * Psi;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScatteredPFTIntegrals(SWGVolume *Oa, int nbfa,
                              SWGVolume *Ob, int nbfb,
                              cdouble Omega, cdouble Q[NUMPFTT])
{
  PFTIData MyData, *Data = &MyData;
  Data->Omega            = real(Omega);
  Data->SameObject       = (Oa==Ob);
  Data->TorqueCenterA    = Oa->Origin;
  Data->TorqueCenterB    = Ob->Origin;

  double Error[2*(NUMPFTT+3)];
  int ncv = CompareBFs(Oa, nbfa, Ob, nbfb);

  double RadiusA=Oa->Faces[nbfa]->Radius;
  double RadiusB=Ob->Faces[nbfb]->Radius;
  double kR = abs(Omega) * fmax(RadiusA, RadiusB);
  bool HighFrequency = ( kR > 5.0 );
  int NumPts = HighFrequency ? ( (ncv > 0) ? 33 : 16 ) : ( (ncv>0) ? 16 : 4 );

  int IDim = 2*(NUMPFTT+3);
  BFBFInt(Oa, nbfa, Ob, nbfb,
          ScatteredPFTIntegrand, (void *)Data, IDim, 
          (double *)Q, Error, NumPts, 0, 0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetEMTPFT(SWGGeometry *G, cdouble Omega, IncField *IF,
                   HVector *JVector, HMatrix *DRMatrix,
                   HMatrix *PFTMatrix, bool Itemize)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NO = G->NumObjects;
  if (    (PFTMatrix==0)
       || (PFTMatrix->NR != NO)
       || (PFTMatrix->NC != NUMPFT)
     )
   ErrExit("invalid PFTMatrix in GetEMTPFT");

  /***************************************************************/
  /* ScatteredPFTT[no] = contributions of object #no to PFTT      */
  /***************************************************************/
  static int NOSave=0;
  static HMatrix **ScatteredPFTT=0, *ExtinctionPFTT=0;
  if (NOSave!=NO)
   { if (ScatteredPFTT)
      { for(int no=0; no<NOSave; no++)
         if (ScatteredPFTT[no]) 
          delete ScatteredPFTT[no];
        free(ScatteredPFTT);
       if (ExtinctionPFTT)
        delete ExtinctionPFTT;
      };
     NOSave=NO;
     ScatteredPFTT=(HMatrix **)mallocEC(NO*sizeof(HMatrix));
     for(int no=0; no<NO; no++)
      ScatteredPFTT[no]=new HMatrix(NO, NUMPFTT);
     ExtinctionPFTT=new HMatrix(NO, NUMPFTT);
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
  memset(DeltaPFTT, 0, NTNO2NQ*sizeof(double));

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
  int TotalBFs    = G->TotalBFs;
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

      cdouble Q[NUMPFTT+3];
      GetScatteredPFTIntegrals(OA, nfa, OB, nfb, Omega, Q);
   
      cdouble u0JJ = ZVAC*GetJJ(JVector, DRMatrix, nbfa, nbfb);
      if (u0JJ==0.0) continue;

      double dPFTT[NUMPFTT];
      if (noa==nob)
       { 
         dPFTT[PFT_PSCAT] = real(u0JJ)*imag(Q[PFT_PSCAT]);
         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] = imag(u0JJ)*imag(Q[nq]);
       }
      else
       {
         dPFTT[PFT_PSCAT] = -1.0*real(u0JJ * II*Q[PFT_PSCAT]);

         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] = -1.0*imag(u0JJ * II * Q[nq]);
       };

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      int Offset = nt*NO2NQ + noa*NONQ + nob*NQ;
      VecPlusEquals(DeltaPFTT + Offset, 0.5, dPFTT, NUMPFTT);

/***************************************************************/
/***************************************************************/
/***************************************************************/
if (UseSymmetry && nbfa!=nbfb)
{
  Q[PFT_XFORCE]*=-1.0;
  Q[PFT_YFORCE]*=-1.0;
  Q[PFT_ZFORCE]*=-1.0;
  Q[PFT_XTORQUE1]*=-1.0;
  Q[PFT_YTORQUE1]*=-1.0;
  Q[PFT_ZTORQUE1]*=-1.0;
  Q[PFT_XTORQUE2]=Q[PFT_XTORQUE2+3];
  Q[PFT_YTORQUE2]=Q[PFT_YTORQUE2+3];
  Q[PFT_ZTORQUE2]=Q[PFT_ZTORQUE2+3];

      cdouble u0JJba = ZVAC*GetJJ(JVector, DRMatrix, nbfb, nbfa);
      if (noa==nob)
       { 
         dPFTT[PFT_PSCAT] = real(u0JJba)*imag(Q[PFT_PSCAT]);
         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] += imag(u0JJba)*imag(Q[nq]);
       }
      else
       {
         dPFTT[PFT_PSCAT] = -1.0*real(u0JJba * II*Q[PFT_PSCAT]);

         for(int nq=PFT_XFORCE; nq<NUMPFTT; nq++)
          dPFTT[nq] = -1.0*imag(u0JJba * II * Q[nq]);
       };

      Offset = nt*NO2NQ + nob*NONQ + noa*NQ;
      VecPlusEquals(DeltaPFTT + Offset, 0.5, dPFTT, NUMPFTT);
};
/***************************************************************/
/***************************************************************/
/***************************************************************/

    }; // end of multithreaded loop

  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  for(int no=0; no<NO; no++)
   ScatteredPFTT[no]->Zero();
  for(int noa=0; noa<NO; noa++)
   for(int nob=0; nob<NO; nob++)
    for(int nq=0; nq<NQ; nq++)
     for(int nt=0; nt<NT; nt++)
      ScatteredPFTT[nob]->AddEntry(noa, nq, DeltaPFTT[ nt*NO2NQ + noa*NONQ + nob*NQ + nq ]);

  /***************************************************************/
  /* get incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   GetExtinctionPFTT(G, JVector, IF, Omega, ExtinctionPFTT);
  else 
   ExtinctionPFTT->Zero();
   
  /***************************************************************/
  /* sum scattered contributions from all objects plus           */
  /* extinction contributions to get total PFT.                  */
  /* note that the formula for total PFT is                      */
  /* Q^{full} = Q^{extinction} - Q^{scattering}                  */
  /* so the scattering contributions enter with a minus sign     */
  /* (except for the scattered power).                           */
  /***************************************************************/
  PFTMatrix->Zero();
  for(int noa=0; noa<NO; noa++)
   for(int nq=0; nq<NUMPFT; nq++)
    { 
      PFTMatrix->AddEntry(noa,nq,ExtinctionPFTT->GetEntry(noa,nq));

      for(int nob=0; nob<NO; nob++)
       { 
         if (nq==PFT_PABS)
          PFTMatrix->AddEntry(noa,nq,-1.0*ScatteredPFTT[nob]->GetEntry(noa,PFT_PSCAT));
         else if (nq==PFT_PSCAT)
          PFTMatrix->AddEntry(noa,nq,+1.0*ScatteredPFTT[nob]->GetEntry(noa,PFT_PSCAT));
         else // force or torque 
          PFTMatrix->AddEntry(noa,nq,-1.0*ScatteredPFTT[nob]->GetEntry(noa,nq));
       };

      if (PFT_XTORQUE<=nq && nq<=PFT_ZTORQUE)
       { PFTMatrix->AddEntry(noa,nq,ExtinctionPFTT->GetEntry(noa,nq+3));
         for(int nob=0; nob<NO; nob++)
          PFTMatrix->AddEntry(noa, nq, -1.0*ScatteredPFTT[nob]->GetEntry(noa,nq+3));
       };
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *ss=getenv("BUFF_ITEMIZE_PFT");
  if (ss && ss[0]=='1')
   Itemize=true;
  if (Itemize)
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
            fprintf(f,"# %i-%i PAbs, PScat, Fxyz, T1xyz T2xyz (object %s)\n",nc,nc+NUMPFTT-1,G->Objects[nob]->Label);
         };
        fprintf(f,"%e %s ",real(Omega),G->Objects[noa]->Label);
        for(int nq=0; nq<NUMPFT; nq++)
         fprintf(f,"%e ",PFTMatrix->GetEntryD(noa,nq));
        for(int nq=0; nq<NUMPFTT; nq++)
         fprintf(f,"%e ",ExtinctionPFTT->GetEntryD(noa,nq));
        for(int nob=0; nob<NO; nob++)
         for(int nq=0; nq<NUMPFTT; nq++)
          fprintf(f,"%e ",ScatteredPFTT[nob]->GetEntryD(noa,nq));
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
