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
        16, 0, 0);
}

/***************************************************************/
/* Return n1, n2 such that nPair = n1*N1 + n1*(n1-1)/2 + (n2-n1)*/
/***************************************************************/
void GetPairIndices(int nPair, int N, int *pn1, int *pn2)
{
  int n1=nPair/N;
  int nDelta = nPair- n1*N - n1*(n1-1)/2;
  while (nDelta<0)
   { n1--;
     nDelta = nPair- n1*N - n1*(n1-1)/2;
   };
  *pn1=n1; 
  *pn2=n1+nDelta;
}

/***************************************************************/
/* PFT[no][nq] = nqth PFT quantity for noth object             */
/***************************************************************/
HMatrix *SWGGeometry::GetPFT(IncField *IF, HVector *JVector,
                             cdouble Omega, HMatrix *PFTMatrix,
                             bool *NeedQuantity, void *opTable)
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

  cdouble IZ = II*ZVAC;
  cdouble IKZ = Omega*IZ;

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
   
      int NPairs    = (OA==OB) ? (NBFA*(NBFA+1)/2) : (NBFA*NBFB);

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
         if (OA==OB)
          GetPairIndices(nPair, NBFA, &nbfA, &nbfB);
         else
          { nbfA      = nPair / NBFB;
            nbfB      = nPair % NBFB;
          };

         if(OA==OB && nbfB==nbfA)
          LogPercent(nPair,NPairs,10);
         else if (OA!=OB && nbfB==0)
          LogPercent(nbfA,NBFA,10);
   
         cdouble G, dG[6];
         G=GetGMatrixElement(OA, nbfA, OB, nbfB, Omega, opTable, dG, NeedQuantity);
         cdouble JJ = conj ( JVector->GetEntry(OffsetA + nbfA) )
                          *( JVector->GetEntry(OffsetB + nbfB) );
   
         double Factor = 0.5;
         if (noA==noB && nbfB > nbfA)
          Factor *= 2.0;
   
         P  += Factor * real ( JJ * IKZ * G    );
         Fx += Factor * imag ( JJ * IZ * dG[0] );
         Fy += Factor * imag ( JJ * IZ * dG[1] );
         Fz += Factor * imag ( JJ * IZ * dG[2] );
         Tx += Factor * imag ( JJ * IZ * dG[3] );
         Ty += Factor * imag ( JJ * IZ * dG[4] );
         Tz += Factor * imag ( JJ * IZ * dG[5] );

       };  // end of multithreaded loop
    
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
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  for(int no=0; no<NumObjects; no++)
   { 
     SWGVolume *O = Objects[no];
     int Offset   = BFIndexOffset[no];
     int NBF      = O->NumInteriorFaces;

     double P2=0.0;
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

           P2 += 0.5*real( JJ * IKZVInv );
         }; // for(int nnz=0; nnz<NNZ; nnz++)

      }; // for(int nbfA=0; nbfA<NBF; nbfA++)
     PFTMatrix->SetEntry(no, 4, P2);

   }; // for(int no=0; no<NumObjects; no++)
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
LogTaskTiming("PFT");
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return PFTMatrix;

}
  
} // namespace buff
