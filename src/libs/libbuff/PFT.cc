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

// numbers of PFT quantities and integrals
#define NUMPFTQS 7
#define NUMPFTIS 7

namespace buff {

typedef struct PFTIntegrandData
 {
   cdouble k;
   bool NeedDerivatives;

 } PFTIntegrandData;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
int MaxTDEvals=100000;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFBF(double *xA, double *bA, double DivbA,
                       double *xB, double *bB, double DivbB,
                       void *UserData, double *I)
{
  PFTIntegrandData Data = (PFTIntegrandData *)UserData;
  cdouble k             = Data->k;
  bool NeedDerivatives  = Data->NeedDerivatives;

  double R[3];
  R[0]=xA[0]-xB[0];
  R[1]=xA[1]-xB[1];
  R[2]=xA[2]-xB[2];
  double r=sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  if ( fabs(r)<1.0e-12 )
   { memset(I, 0, 2*NUMPFTIS*sizeof(double));
     return;
   };

  cdouble ScalarFactor = VecDot(bA, bB) + DivbA*DivbB/k2;
  cdouble ExpFac=exp(II*k*r)/(4.0*M_PI*r);
  cdouble *zI = (cdouble *)I;
  zI[0] = ScalarFactor * ExpFac;
  if (NeedDerivatives)
   { cdouble Psi=(II*k*r-1.0)*ExpFac/(r*r);
     zI[1] = ScalarFactor * R[0] * Psi;
     zI[2] = ScalarFactor * R[1] * Psi;
     zI[3] = ScalarFactor * R[2] * Psi;
     zI[4] = 0.0;
     zI[5] = 0.0;
     zI[6] = 0.0;
   };

}

/***************************************************************/
/* compute PFT integrals between pairs of SWG basis functions  */
/***************************************************************/
void GetPFTIntegrals_BFBF(SWGVolume *OA, int nbfA,
                          SWGVolume *OB, int nbfB,
                          cdouble Omega, cdouble IPFT[NUMPFTIS])
{ 
  double rRel=0.0;
  int ncv= (OA==OB) ? CompareBFs(OA, nbfA, OB, nbfB, &rRel) : 0;
  
  PFTIntegrandData MyData, *Data=&MyData;
  Data->k               = Omega;
  Data->NeedDerivatives = true;

  if (ncv==0)
   { 
     /* if there are no common vertices, use low-order cubature */
     cdouble Error[NUMPFTIS];
     BFBFInt(OA, nbfA, OB, nbfB, PFTIntegrand_BFBF, (void *)Data,
             2*NUMPFTIS, (double *)IPFT, (double *)Error, 4, 0, 0);
   }
  else 
   { 
     /* use Taylor-Duffy for pairs of tetrahedra with 3 or 4 */
     /* common vertices, and high-order cubature otherwise   */
     SWGFace *FA = OA->Faces[nbfA];
     SWGFace *FB = OB->Faces[nbfB];
     memset(IPFT, 0, NUMPFTIS*sizeof(cdouble));
     for(int SignA=1; SignA>=-1; SignA-=2)
      for(int SignB=1; SignB>=-1; SignB-=2)
       { 
         int ntA = (SignA==1) ? FA->iPTet  : FA->iMTet;
         int iQA = (SignA==1) ? FA->PIndex : FA->MIndex;
         int ntB = (SignB==1) ? FB->iPTet  : FB->iMTet;
         int iQB = (SignB==1) ? FB->PIndex : FB->MIndex;

         int ncv=CompareTets(OA, ntA, OB, ntB);
         cdouble Result[NUMPFTIS], Error[NUMPFTIS];
    
         if (ncv<=3)
          TetTetInt(OA, ntA, iQA, 1.0, OB, ntB, iQB, 1.0,
                    PFTIntegrand_BFBF, (void *)Data, 2*NUMPFTIS, 
                    (double *)Result, (double *)Error, 
                    33, 0, 0);
         else
          TetTetInt_TD(OA, ntA, iQA, OB, ntB, iQB,
                       PFTIntegrand_BFBF, (void *)Data, 2*NUMPFTIS,
                       (double *)Result, (double *)Error,
                       0, MaxTDEvals, 1.0e-8);

         for(int n=0; n<NUMPFTIS; n++)
          IPFT[n] += SignA*SignB*Result[n];
       };
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PFTIntegrand_BFInc(double *x, double *b, double Divb,
                        void *UserData, double *I)
{
  PFTIntegrandData *PFTIData=(PFTIntegrandData *)UserData;

  cdouble k    = PFTIData->k;
  IncField *IF = PFTIData->IF; 

  // get fields and derivatives at eval point by finite-differencing
  cdouble EH[6], EHP[6], EHM[6], dEH[3][6];
  IF->GetFields(x, EH);
  for(int Mu=0; Mu<3; Mu++)
   { 
     double Delta = 1.0e-4 / abs(k);

     x[Mu] += Delta;
     IF->GetFields(x, EHP);
     x[Mu] -= 2.0*Delta;
     IF->GetFields(x, EHM);
     x[Mu] += Delta;

     for(int Nu=0; Nu<6; Nu++)
      dEH[Mu][Nu] = (EHP[Nu]-EHM[Nu])/(2.0*Delta);
   };

  cdouble *zI = (cdouble *)I;
  memset(zI, 0, NUMPFTIS*sizeof(cdouble));
  for(Mu=0; Mu<3; Mu++)
   { zI[0] += b[Mu]*EH[Mu];
     zI[1] += b[Mu]*dEH[0][Mu];
     zI[2] += b[Mu]*dEH[1][Mu];
     zI[3] += b[Mu]*dEH[2][Mu];
     zI[4] += 0.0;
     zI[5] += 0.0;
     zI[6] += 0.0;
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

  cdouble Error[NUMPFTIS];
  BFInt(O, nbf, PFTIntegrand_BFInc, (void *)PFTIData,
        2*NUMPFTIS, (double *)IPFT, (double *)Error,
        16, 0, 0);
}

/***************************************************************/
/* PFT[no][nq] = nqth PFT quantity for noth object             */
/***************************************************************/
HMatrix *SWGGeometry::GetPFT(IncField *IF, HVector *JVector,
                             cdouble Omega, HMatrix *PFTMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( (char *s=getenv("BUFF_MAXTDEVALS")) )
   { sscanf(s,"%i",&MaxTDEvals);
     Log("Setting MaxTDEvals=%i.",MaxTDEvals);
   };

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
  int NumThreads = GetNumThreads();
  for(int noA=0; noA<NumObjects; noA++)
   for(int noB=noA; noB<NumObjects; noB++)
    { 
      SWGObject *OA = Objects[noA];
      int OffsetA   = BFIndexOffset[noA];
      int NBFA      = OA->NumInteriorFaces;

      SWGObject *OB = Objects[noB];
      int OffsetB   = BFIndexOffset[noB];
      int NBFB      = OB->NumInteriorFaces;

      int NBFP;
      bool UseSymmetry = false;
      if (UseSymmetry)
       NBFP = OA==OB ? (NBFA*(NBFA+1)/2) : (NBFA*NBFB);
      else 
       NBFP = NBFA*NBFB;

      /*--------------------------------------------------------------*/
      /*- multithreaded loop over basis functions on OA, OB-----------*/
      /*--------------------------------------------------------------*/
      double P=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Tx=0.0, Ty=0.0, Tz=0.0; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:P, Fx, Fy, Fz, Tx, Ty, Tz)
#endif
      for(int nbfp=0; nbfp<NBFP; nbfp++)
       { 
         int nbfA, nbfB;
         if (UseSymmetry)
          { nbfA      = NBFP / NBFB;
            nbfB      = NBFP % NBFB;
          }
         else
          { nbfA      = NBFP / NBFB;
            nbfB      = NBFP % NBFB;
          };

         double IPFT[7];
         GetPFTIntegrals(OA, nbfA, OB, nbfB, Omega, IPFT);
         cdouble JJ = conj ( JVector->GetEntry(OffsetA + nbfA) )
                          *( JVector->GetEntry(OffsetB + nbfB) );

         P  += real( JJ * IPFT[0] );
         Fx += imag( JJ * IPFT[1] );
         Fy += imag( JJ * IPFT[2] );
         Fz += imag( JJ * IPFT[3] );
         Tx += imag( JJ * IPFT[4] );
         Ty += imag( JJ * IPFT[5] );
         Tz += imag( JJ * IPFT[6] );
       };

      /*--------------------------------------------------------------*/
      /*- accumulate PFT contributions for this pair of objects       */
      /*--------------------------------------------------------------*/
      double Factor = 0.5;
      if (noA==noB && nbfB > nbfA)
       Factor = 1.0;

      PFTMatrix->AddEntry(noA, 0, Factor * P  );
      if (noB!=noA)
       PFTMatrix->AddEntry(noB, 0, Factor * P );

      Factor/ = real(Omega);
      PFTMatrix->AddEntry(noA, 1, Factor * Fx );
      PFTMatrix->AddEntry(noA, 2, Factor * Fy );
      PFTMatrix->AddEntry(noA, 3, Factor * Fz );
      PFTMatrix->AddEntry(noA, 4, Factor * Tx );
      PFTMatrix->AddEntry(noA, 5, Factor * Ty );
      PFTMatrix->AddEntry(noA, 6, Factor * Tz );
      if (noB>noA)
       { PFTMatrix->AddEntry(noB, 1, -Factor * Fx );
         PFTMatrix->AddEntry(noB, 2, -Factor * Fy );
         PFTMatrix->AddEntry(noB, 3, -Factor * Fz );
         PFTMatrix->AddEntry(noB, 4, -Factor * Tx );
         PFTMatrix->AddEntry(noB, 5, -Factor * Ty );
         PFTMatrix->AddEntry(noB, 6, -Factor * Tz );
       };

    };  // for(int noA=0:NumObjects, noB=noA:NumObjects

  /***************************************************************/
  /* add incident-field contributions ****************************/
  /***************************************************************/
  if (IF)
   { 
      for(int no=0; no<NO; no++)
       { 
         SWGObject *O = Objects[no];
         int Offset   = BFIndexOffset[no];
         int NBF      = O->NumInteriorFaces;

         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         /*--------------------------------------------------------------*/
         double P=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Tx=0.0, Ty=0.0, Tz=0.0; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:P, Fx, Fy, Fz, Tx, Ty, Tz)
#endif
         for(int nbf=0; nbf<O->NumInteriorFaces; nbf++)
          { 
            double IPFT[7];
            GetPFTIntegrals(O, nbf, IF, Omega, IPFT);
            cdouble J = JVector->GetEntry(Offset + nbf);
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
         Factor/=Omega;
         PFTMatrix->AddEntry(no, 1, Factor * Fx );
         PFTMatrix->AddEntry(no, 2, Factor * Fy );
         PFTMatrix->AddEntry(no, 3, Factor * Fz );
         PFTMatrix->AddEntry(no, 4, Factor * Tx );
         PFTMatrix->AddEntry(no, 5, Factor * Ty );
         PFTMatrix->AddEntry(no, 6, Factor * Tz );
       };

   }; // if (IF)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return PFTMatrix;

}
  
} // namespace buff
