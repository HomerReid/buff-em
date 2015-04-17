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

// numbers of PFT quantities and integrals
#define NUMPFTQS 7
#define NUMPFTIS 7

int GetVInvAndImEpsEntries(SWGVolume *V, int nfA,
                           cdouble Omega, int Indices[7],
                           cdouble VInvEntries[7],
                           double ImEpsEntries[7]);

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
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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
     double Delta = 1.0e-4 / abs(k);

     xTweaked[Mu] += Delta;
     IF->GetFields(xTweaked, EHP);
     xTweaked[Mu] -= 2.0*Delta;
     IF->GetFields(xTweaked, EHM);
     xTweaked[Mu] += Delta;

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
/* Return n1, n2 such that nPair = n1*N1 - n1*(n1-1)/2 + (n2-n1)*/
/***************************************************************/
void GetPairIndices(int nPair, int N, int *pn1, int *pn2)
{
  int n1=nPair/N;
  int nDelta = nPair - n1*N + n1*(n1-1)/2;
  while ( (n1+nDelta) >= N )
   { n1++;
     nDelta = nPair - n1*N + n1*(n1-1)/2;
   };
  *pn1=n1; 
  *pn2=n1+nDelta;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SWGGeometry::GetSparsePFT(HVector *JVector,
                                   cdouble Omega, HMatrix *PFTMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMatrix!=0 && (PFTMatrix->NR!=NumObjects || PFTMatrix->NC!=NUMPFTQS) )
   { delete PFTMatrix;
     PFTMatrix=0;
   };
  if (PFTMatrix==0)
   PFTMatrix= new HMatrix(NumObjects, NUMPFTQS);

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
        int nbfBList[7];
        cdouble VInvList[7];
        double ImEpsList[7];
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

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
        cdouble BxImVB[7][3];
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
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

      }; // for(int nbfA=0; nbfA<NBF; nbfA++)
     PFTMatrix->SetEntry(no, 0, P2);
     PFTMatrix->SetEntry(no, 1, T[0]);
     PFTMatrix->SetEntry(no, 2, T[1]);
     PFTMatrix->SetEntry(no, 3, T[2]);

   }; // for(int no=0; no<NumObjects; no++)

  /***************************************************************/
  /* convert force/torque values to units of nanoNewtons         */
  /***************************************************************/
  #define TENTHIRDS (10.0/3.0)
  for(int no=0; no<NumObjects; no++)
   for(int Mu=1; Mu<=6; Mu++)
    PFTMatrix->SetEntry(no,Mu, TENTHIRDS*PFTMatrix->GetEntry(no,Mu) );

  return PFTMatrix;

}

/***************************************************************/
/* PFT[no][nq] = nqth PFT quantity for noth object             */
/***************************************************************/
HMatrix *SWGGeometry::GetDensePFT(IncField *IF, HVector *JVector,
                                  cdouble Omega, HMatrix *PFTMatrix,
                                  bool *NeedQuantity)
{ 
  bool DefaultNeedQuantity[6]={true, true, true, true, true, true};
  if (NeedQuantity==0) NeedQuantity=DefaultNeedQuantity;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMatrix!=0 && (PFTMatrix->NR!=NumObjects || PFTMatrix->NC!=NUMPFTQS) )
   { delete PFTMatrix;
     PFTMatrix=0;
   };
  if (PFTMatrix==0)
   PFTMatrix= new HMatrix(NumObjects, NUMPFTQS);
  PFTMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool UseSymmetry=false;
  char *s=getenv("BUFF_PFT_SYMMETRY");
  if (s && s[0]=='1')
   UseSymmetry=true;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int noA=0; noA<NumObjects; noA++)
   for(int noB=noA; noB<NumObjects; noB++)
    { 
      SWGVolume *OA = Objects[noA];
      int OffsetA   = BFIndexOffset[noA];
      int NBFA      = OA->NumInteriorFaces;
   
      SWGVolume *OB = Objects[noB];
      int OffsetB   = BFIndexOffset[noB];
      int NBFB      = OB->NumInteriorFaces;

      int NPairs    = (UseSymmetry && OA==OB) ? (NBFA*(NBFA+1)/2) : (NBFA*NBFB);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
const char *TaskNames[]={ "NCV0", "NCV1", "NCV2", "NCV3", "NCV4","IF  ",0};
InitTaskTiming( TaskNames );
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   
      /*--------------------------------------------------------------*/
      /*- multithreaded loop over basis functions on OA, OB-----------*/
      /*--------------------------------------------------------------*/
      double P=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Tx=0.0, Ty=0.0, Tz=0.0;
      Log("Computing PFT (%i,%i)...",noA,noB);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
feenableexcept(FE_INVALID | FE_OVERFLOW);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#ifdef USE_OPENMP
int NumThreads = GetNumThreads();
Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:P, Fx, Fy, Fz, Tx, Ty, Tz)
#endif
      for(int nPair=0; nPair<NPairs; nPair++)
       { 
         int nbfA=0, nbfB=0;
         if (UseSymmetry && OA==OB)
          GetPairIndices(nPair, NBFA, &nbfA, &nbfB);
         else
          { nbfA      = nPair / NBFB;
            nbfB      = nPair % NBFB;
          };

         if(UseSymmetry && OA==OB && nbfB==nbfA)
          LogPercent(nPair,NPairs,10);
         else if ( (!UseSymmetry || OA!=OB) && nbfB==0)
          LogPercent(nbfA,NBFA,10);

         cdouble JJ = conj ( JVector->GetEntry(OffsetA + nbfA) )
                          *( JVector->GetEntry(OffsetB + nbfB) );
         if (JJ==0.0) continue;
   
         cdouble G, dG[6];
         G=GetGMatrixElement(OA, nbfA, OB, nbfB, Omega, Cache, dG, NeedQuantity);
   
         if (UseSymmetry)
          { 
            double Factor = ( OA==OB && nbfB == nbfA ) ? 0.5 : 1.0;
            P  -= Factor * real ( JJ ) * imag(G);
            Fx -= imag(JJ) * imag(dG[0]);
            Fy -= imag(JJ) * imag(dG[1]);
            Fz -= imag(JJ) * imag(dG[2]);
            Tx -= imag(JJ) * imag(dG[3]);
            Ty -= imag(JJ) * imag(dG[4]);
            Tz -= imag(JJ) * imag(dG[5]);
          }
         else
          { P  += 0.5 * real ( II*JJ*G     );
            Fx += 0.5 * imag ( II*JJ*dG[0] );
            Fy += 0.5 * imag ( II*JJ*dG[1] );
            Fz += 0.5 * imag ( II*JJ*dG[2] );
            Tx += 0.5 * imag ( II*JJ*dG[3] );
            Ty += 0.5 * imag ( II*JJ*dG[4] );
            Tz += 0.5 * imag ( II*JJ*dG[5] );
          };

       };  // end of multithreaded loop
      P  *= ZVAC*real(Omega);
      Fx *= ZVAC;
      Fy *= ZVAC;
      Fz *= ZVAC;
      Tx *= ZVAC;
      Ty *= ZVAC;
      Tz *= ZVAC;
    
      /*--------------------------------------------------------------*/
      /*- accumulate PFT contributions for this pair of objects       */
      /*--------------------------------------------------------------*/
      PFTMatrix->AddEntry(noA, 0, P  );
      PFTMatrix->AddEntry(noA, 1, Fx );
      PFTMatrix->AddEntry(noA, 2, Fy );
      PFTMatrix->AddEntry(noA, 3, Fz );
      PFTMatrix->AddEntry(noA, 4, Tx );
      PFTMatrix->AddEntry(noA, 5, Ty );
      PFTMatrix->AddEntry(noA, 6, Tz );
      if (noB>noA)
       { PFTMatrix->AddEntry(noB, 0, P );
         PFTMatrix->AddEntry(noB, 1, -1.0 * Fx );
         PFTMatrix->AddEntry(noB, 2, -1.0 * Fy );
         PFTMatrix->AddEntry(noB, 3, -1.0 * Fz );
         PFTMatrix->AddEntry(noB, 4, -1.0 * Tx );
         PFTMatrix->AddEntry(noB, 5, -1.0 * Ty );
         PFTMatrix->AddEntry(noB, 6, -1.0 * Tz );
       };
   
   };  // for(int noA=0:NumObjects, noB=noA:NumObjects

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
int NumThreads = GetNumThreads();
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
         Factor/=real(Omega);
         PFTMatrix->AddEntry(no, 1, Factor * Fx );
         PFTMatrix->AddEntry(no, 2, Factor * Fy );
         PFTMatrix->AddEntry(no, 3, Factor * Fz );
         PFTMatrix->AddEntry(no, 4, Factor * Tx );
         PFTMatrix->AddEntry(no, 5, Factor * Ty );
         PFTMatrix->AddEntry(no, 6, Factor * Tz );
       };

   }; // if (IF)
  Elapsed = Secs() - Elapsed;
  AddTaskTiming(5,Elapsed);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
LogTaskTiming("DensePFT");
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /***************************************************************/
  /* convert force/torque values to units of nanoNewtons         */
  /***************************************************************/
  #define TENTHIRDS (10.0/3.0)
  for(int no=0; no<NumObjects; no++)
   for(int Mu=1; Mu<=6; Mu++)
    PFTMatrix->SetEntry(no,Mu, TENTHIRDS*PFTMatrix->GetEntry(no,Mu) );

  return PFTMatrix;

}
  
} // namespace buff
