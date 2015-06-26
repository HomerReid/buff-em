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
 * OPFT.cc    -- libbuff class methods for computing power, force,
 *            -- and torque in classical deterministic scattering
 *            -- problems using the "overlap" formalism
 
 * homer reid -- 1/2015
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
namespace buff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
cdouble GetJJ(HVector *JVector, HMatrix *Rytov, int nbfa, int nbfb);
SWGVolume *ResolveNBF(SWGGeometry *G, int nbf, int *pno, int *pnf);
void Invert3x3Matrix(cdouble M[3][3], cdouble W[3][3]);
int GetVInvEntries(SWGVolume *V, int nfA, cdouble Omega,
                   int Indices[MAXOVERLAP], 
                   cdouble VInvEntries[MAXOVERLAP]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
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
HMatrix *GetOPFT(SWGGeometry *G, cdouble Omega,
                 HVector *JVector, HMatrix *Rytov,
                 HMatrix *PFTMatrix)
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
  for(int no=0; no<G->NumObjects; no++)
   { 
     SWGVolume *O = G->Objects[no];
     int Offset   = G->BFIndexOffset[no];
     int NBF      = O->NumInteriorFaces;

     double P=0.0;
     double T[3]={0.0, 0.0, 0.0};
     for(int nbfA=0; nbfA<NBF; nbfA++)
      { 
        int nbfBList[MAXOVERLAP];
        cdouble VInvList[MAXOVERLAP];
        cdouble BxImVB[MAXOVERLAP][3];
        int NNZ=GetVInvEntries(O, nbfA, Omega, nbfBList, VInvList);
        int NNZ2=GetBxImVB(O, nbfA, Omega, nbfBList, BxImVB);
        if (NNZ2!=NNZ) 
         ErrExit("%s:%i: internal error",__FILE__,__LINE__);

        for(int nnz=0; nnz<NNZ; nnz++)
         { 
           int nbfB = nbfBList[nnz];
           cdouble JJ=GetJJ(JVector, Rytov, Offset+nbfA, Offset+nbfB);
           cdouble IKZVInv = II*Omega*ZVAC*VInvList[nnz];
           P -= 0.5*real( JJ * IKZVInv );
           T[0] += 0.5*TENTHIRDS*ZVAC*real(JJ*BxImVB[nnz][0]);
           T[1] += 0.5*TENTHIRDS*ZVAC*real(JJ*BxImVB[nnz][1]);
           T[2] += 0.5*TENTHIRDS*ZVAC*real(JJ*BxImVB[nnz][2]);
         }

      }; // for(int nbfA=0; nbfA<NBF; nbfA++)

     PFTMatrix->SetEntry(no, PFT_PABS,   P);
     PFTMatrix->SetEntry(no, PFT_XFORCE, T[0]);
     PFTMatrix->SetEntry(no, PFT_YFORCE, T[1]);
     PFTMatrix->SetEntry(no, PFT_ZFORCE, T[2]);

   }; // for(int no=0; no<NumObjects; no++)

  return PFTMatrix;
}
  
} // namespace buff
