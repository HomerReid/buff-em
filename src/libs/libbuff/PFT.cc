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
 * PFT.cc        -- libbuff class methods for computing power, force,
 *               -- and torque in classical deterministic scattering
 *               -- problems
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

int GetVInvAndImEpsEntries(SWGVolume *V, int nfA,
                           cdouble Omega, int Indices[MAXOVERLAP],
                           cdouble VInvEntries[MAXOVERLAP],
                           double ImEpsEntries[MAXOVERLAP]);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void Invert3x3Matrix(cdouble M[3][3], cdouble W[3][3]);

typedef struct VII2Data
 { 
   double *QA;
   double PreFacA;
   double *QB;
   double PreFacB;
   cdouble Omega;
   IHAIMatProp *MP;

 } VII2Data;

void VIntegrand2(double *x, double *b, double DivB, 
                  void *UserData, double *I)
{
  (void) DivB;
  (void) b;

  VII2Data *Data   = (VII2Data *)UserData;
  double *QA      = Data->QA;
  double PreFacA  = Data->PreFacA;
  double *QB      = Data->QB;
  double PreFacB  = Data->PreFacB;
  cdouble Omega   = Data->Omega;
  IHAIMatProp *MP = Data->MP;

  cdouble Eps[3][3], Y[3][3];
  MP->GetEps( Omega, x, Eps );

  Eps[0][0] -= 1.0;
  Eps[1][1] -= 1.0;
  Eps[2][2] -= 1.0;
  Invert3x3Matrix(Eps, Y);

  double FA[3], FB[3];
  FA[0] = PreFacA * (x[0] - QA[0]);
  FA[1] = PreFacA * (x[1] - QA[1]);
  FA[2] = PreFacA * (x[2] - QA[2]);
  FB[0] = PreFacB * (x[0] - QB[0]);
  FB[1] = PreFacB * (x[1] - QB[1]);
  FB[2] = PreFacB * (x[2] - QB[2]);

  cdouble YFB[3];
  cdouble Omega2=Omega*Omega;
  YFB[0] = -(Y[0][0]*FB[0] + Y[0][1]*FB[1] + Y[0][2]*FB[2]) / Omega2;
  YFB[1] = -(Y[1][0]*FB[0] + Y[1][1]*FB[1] + Y[1][2]*FB[2]) / Omega2;
  YFB[2] = -(Y[2][0]*FB[0] + Y[2][1]*FB[1] + Y[2][2]*FB[2]) / Omega2;

  cdouble *zI=(cdouble *)I;
  zI[0] = FA[1]*YFB[2] - FA[2]*YFB[1];
  zI[1] = FA[2]*YFB[0] - FA[0]*YFB[2];
  zI[2] = FA[0]*YFB[1] - FA[1]*YFB[0];
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetBxImVB(SWGVolume *V, int nfA, cdouble Omega,
              int Indices[7], cdouble BxImVB[7][3])
{
  Indices[0]=nfA;
  int NNZ=1;

  SWGFace *FA = V->Faces[nfA];
  VII2Data MyVII2Data, *Data=&MyVII2Data;
  Data->Omega = Omega;
  Data->MP    = V->MP;

  /*--------------------------------------------------------------*/
  /* handle interactions between bfs #nfA and #nfB, where         */
  /* nfB runs over the 4 faces of the positive tetrahedron of nfA */
  /*--------------------------------------------------------------*/
  int nt        = FA->iPTet;
  SWGTet *T     = V->Tets[nt];
  Data->QA      = V->Vertices + 3*(FA->iQP);
  Data->PreFacA = FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB >= V->NumInteriorFaces) continue;

     SWGFace *FB = V->Faces[nfB];
     if ( FB->iPTet == nt )
      { Data->QB      = V->Vertices + 3*(FB->iQP);
        Data->PreFacB = FB->Area / (3.0*T->Volume);
      }
     else
      { Data->QB      = V->Vertices + 3*(FB->iQM);
        Data->PreFacB = -1.0*FB->Area / (3.0*T->Volume);
      };

     cdouble I[3], E[3];
     TetInt(V, nt, 0, 1.0, VIntegrand2, (void *)Data,
            6, (double *)I, (double *)E, 33, 0, 1.0e-4);

     if (nfB==nfA)
      { BxImVB[0][0] = I[0];
        BxImVB[0][1] = I[1];
        BxImVB[0][2] = I[2];
      }
     else
      { BxImVB[NNZ][0] = I[0];
        BxImVB[NNZ][1] = I[1];
        BxImVB[NNZ][2] = I[2];
        Indices[NNZ]    = nfB;
        NNZ++;
      };

   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  nt            = FA->iMTet;
  T             = V->Tets[nt];
  Data->QA      = V->Vertices + 3*(FA->iQM);
  Data->PreFacA = -1.0*FA->Area / (3.0*T->Volume);
  for(int ifB=0; ifB<4; ifB++)
   { 
     int nfB = T->FI[ifB];
     if (nfB >= V->NumInteriorFaces) continue;

     SWGFace *FB = V->Faces[nfB];
     if ( FB->iPTet == nt )
      { Data->QB      = V->Vertices + 3*(FB->iQP);
        Data->PreFacB = FB->Area / (3.0*T->Volume);
      }
     else
      { Data->QB      = V->Vertices + 3*(FB->iQM);
        Data->PreFacB = -1.0*FB->Area / (3.0*T->Volume);
      };

     cdouble I[3], E[3];
     TetInt(V, nt, 0, 1.0, VIntegrand2, (void *)Data,
            6, (double *)I, (double *)E, 33, 0, 0);

     if (nfB==nfA)
      { BxImVB[0][0] += I[0];
        BxImVB[0][1] += I[1];
        BxImVB[0][2] += I[2];
      }
     else
      { BxImVB[NNZ][0] = I[0];
        BxImVB[NNZ][1] = I[1];
        BxImVB[NNZ][2] = I[2];
        Indices[NNZ]   = nfB;
        NNZ++;
      };

   }; // for(int ifB=0; ifB<4; ifB++)

  return NNZ;

} // GetBxImVB

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::GetSparsePFT(HVector *JVector,
                                   cdouble Omega, HMatrix *PFTMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMatrix!=0 && (PFTMatrix->NR!=NumObjects || PFTMatrix->NC!=NUMPFT) )
   { delete PFTMatrix;
     PFTMatrix=0;
   };
  if (PFTMatrix==0)
   PFTMatrix= new HMatrix(NumObjects, NUMPFT);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTMatrix->Zero();
  for(int no=0; no<NumObjects; no++)
   { 
     SWGVolume *O = Objects[no];
     int Offset   = BFIndexOffset[no];
     int NBF      = O->NumInteriorFaces;

     double P2=0.0;
     double T[3]={0.0, 0.0, 0.0};
     for(int nbfA=0; nbfA<NBF; nbfA++)
      { 
        int nbfBList[MAXOVERLAP];
        cdouble VInvList[MAXOVERLAP];
        double ImEpsList[MAXOVERLAP];
        int NNZ=GetVInvAndImEpsEntries(O, nbfA, Omega, nbfBList,
                                       VInvList, ImEpsList);

        for(int nnz=0; nnz<NNZ; nnz++)
         { 
           int nbfB = nbfBList[nnz];
           cdouble JJ = conj ( JVector->GetEntry(Offset + nbfA) )
                            *( JVector->GetEntry(Offset + nbfB) );
           cdouble IKZVInv = II*Omega*ZVAC*VInvList[nnz];
           P2 -= 0.5*real( JJ * IKZVInv );
         }; // for(int nnz=0; nnz<NNZ; nnz++)

        cdouble BxImVB[MAXOVERLAP][3];
        NNZ=GetBxImVB(O, nbfA, Omega, nbfBList, BxImVB);
        for(int nnz=0; nnz<NNZ; nnz++)
         { 
           int nbfB = nbfBList[nnz];
           cdouble JJ = conj ( JVector->GetEntry(Offset + nbfA) )
                            *( JVector->GetEntry(Offset + nbfB) );
           T[0] += 0.5*ZVAC*real(JJ*BxImVB[nnz][0]);
           T[1] += 0.5*ZVAC*real(JJ*BxImVB[nnz][1]);
           T[2] += 0.5*ZVAC*real(JJ*BxImVB[nnz][2]);
         }

      }; // for(int nbfA=0; nbfA<NBF; nbfA++)
     PFTMatrix->SetEntry(no, 0, P2);
     PFTMatrix->SetEntry(no, 2, T[0]);
     PFTMatrix->SetEntry(no, 3, T[1]);
     PFTMatrix->SetEntry(no, 4, T[2]);

   }; // for(int no=0; no<NumObjects; no++)

  /***************************************************************/
  /* convert force/torque values to units of nanoNewtons         */
  /***************************************************************/
  #define TENTHIRDS (10.0/3.0)
  for(int no=0; no<NumObjects; no++)
   for(int Mu=2; Mu<=7; Mu++)
    PFTMatrix->SetEntry(no,Mu, TENTHIRDS*PFTMatrix->GetEntry(no,Mu) );

  return PFTMatrix;

}

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
#define TENTHIRDS (10.0/3.0)
HMatrix *SWGGeometry::GetDensePFT(IncField *IF,
                                  HVector *JVector, HVector *RHSVector,
                                  cdouble Omega, HMatrix *PFTMatrix,
                                  bool *NeedFT)
{ 
  bool DefaultNeedFT[6]={true, true, true, true, true, true};
  if (NeedFT==0) NeedFT=DefaultNeedFT;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMatrix!=0 && (PFTMatrix->NR!=NumObjects || PFTMatrix->NC!=NUMPFT) )
   { delete PFTMatrix;
     PFTMatrix=0;
   };
  if (PFTMatrix==0)
   PFTMatrix= new HMatrix(NumObjects, NUMPFT);
  PFTMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (RHSVector)
   { cdouble PreFactor = -II*Omega*ZVAC;
     for(int no=0, nbf=0; no<NumObjects; no++)
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
  int NO = NumObjects;
  int NONQ = NO*NQ;
  double *DeltaPFT = (double *)mallocEC(NumThreads*NONQ*sizeof(double));

  /*--------------------------------------------------------------*/
  /*- multithreaded loop over all basis functions in all volumes -*/
  /*--------------------------------------------------------------*/
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

      int noa=0;
      while ( (noa<NO-1) && nbfa > BFIndexOffset[noa] )
       noa++;
      SWGVolume *OA = Objects[noa];
      int nfa       = nbfa - BFIndexOffset[noa];

      int nob=0;
      while ( (nob<NO-1) && nbfb > BFIndexOffset[nob] )
       nob++;
      SWGVolume *OB = Objects[nob];
      int nfb       = nbfb - BFIndexOffset[nob];
   
      cdouble JJ = conj ( JVector->GetEntry(nbfa) )
                       *( JVector->GetEntry(nbfb) );
      if (JJ==0.0) continue;
   
      cdouble G, dG[6];
      G=GetGMatrixElement(OA, nfa, OB, nfb, Omega, Cache, dG, NeedFT);

      int nt=0;
#ifdef USE_OPENMP
      nt=omp_get_thread_num();
#endif
      
       int Offset = nt*NONQ + noa*NQ;
       DeltaPFT[ Offset + 0 ] += 0.5*real( IKZ*JJ*G );
       for(int Mu=0; Mu<6; Mu++)
        DeltaPFT[ Offset + 2 + Mu ] += 0.5*TENTHIRDS*imag( IKZ*JJ*dG[Mu] ) / real(Omega);

       if (nbfb>nbfa)
        { Offset = nt*NONQ + nob*NQ;
          DeltaPFT[ Offset + 0 ] += 0.5*real( IKZ*conj(JJ*G) );
          for(int Mu=0; Mu<6; Mu++)
           DeltaPFT[ Offset + 2 + Mu ] += 0.5*TENTHIRDS*imag( IKZ*conj(JJ*dG[Mu]) ) / real(Omega);
        };

    }; // end of multithreaded loop
  
  /*--------------------------------------------------------------*/
  /*- accumulate contributions of all threads                     */
  /*--------------------------------------------------------------*/
  for(int no=0; no<NumObjects; no++)
   for(int nq=0; nq<NQ; nq++)
    for(int nt=0; nt<NumThreads; nt++)
     PFTMatrix->AddEntry(no, nq, DeltaPFT[ nt*NONQ + no*NQ + nq ]);

  /***************************************************************/
  /* add incident-field contributions ****************************/
  /***************************************************************/
  double Elapsed=Secs();
  if (IF)
   { 
      for(int no=0; no<NumObjects; no++)
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
  AddTaskTiming(5,Elapsed);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int no=0; no<NumObjects; no++)
   PFTMatrix->AddEntry(no, 1, -1.0*PFTMatrix->GetEntry(no,0));

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
LogTaskTiming("DensePFT");
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  return PFTMatrix;

}
  
} // namespace buff
