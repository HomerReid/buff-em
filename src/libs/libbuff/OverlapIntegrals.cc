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
 * OverlapIntegrals.cc -- libbuff routines for computing various
 *                     -- types of overlap integrals between
 *                     -- SWG basis functions
 *
 * homer reid          -- 7/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libSGJC.h"
#include "libTriInt.h"
#include "libscuff.h"
#include "libbuff.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

typedef void (*OverlapIntegrand)(double x[3],
                                 double bA[3], double DivbA,
                                 double bB[3], double DivbB,
                                 SVTensor *SVT, cdouble Omega,
                                 void *UserData, double *Integrand);

namespace buff {

void Invert3x3Matrix(cdouble M[3][3], cdouble W[3][3]);

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetThetaFactor(double Omega, double T);
#if 0
#define BOLTZMANNK 4.36763e-4
double GetThetaFactor(double Omega, double T)
{ 
  if (T==0.0)
   return 0.0;
  return Omega / ( exp( Omega/(BOLTZMANNK*T) ) - 1.0 );
}
#endif

/***************************************************************/
/* fdim = 5 ****************************************************/
/***************************************************************/
void OverlapIntegrand_VVInvSigma(double x[3],
                                 double bA[3],     double DivbAlpha,
                                 double bB[3],     double DivbBeta,
                                 SVTensor *EpsSVT, cdouble Omega,
                                 void *UserData,   double *I)
{
  (void) DivbAlpha; 
  (void) DivbBeta; // unused

  cdouble EpsM1[3][3], InvEpsM1[3][3];
  EpsSVT->Evaluate( Omega, x, EpsM1 );
  EpsM1[0][0] -= 1.0;
  EpsM1[1][1] -= 1.0;
  EpsM1[2][2] -= 1.0;
  Invert3x3Matrix(EpsM1, InvEpsM1);

  double Theta=1.0;
  SVTensor *TemperatureSVT = (SVTensor *)UserData;
  if (TemperatureSVT)
   { 
     cdouble TT[3][3];
     TemperatureSVT->Evaluate(0,x,TT);
     double T=real(TT[0][0]);
     Theta = GetThetaFactor(real(Omega), T) - GetThetaFactor(real(Omega), 0.0);
   };

  cdouble V=0.0, VInv=0.0;
  double Sigma=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { V     += bA[Mu]*EpsM1[Mu][Nu]*bB[Nu];
      VInv  += bA[Mu]*InvEpsM1[Mu][Nu]*bB[Nu];
      Sigma += Theta*bA[Mu]*imag(EpsM1[Mu][Nu])*bB[Nu];
    };
  V     *= -1.0*Omega*Omega;
  VInv  *= -1.0 / (Omega*Omega);
  Sigma *= 2.0*real(Omega) / (M_PI*ZVAC);
 
  I[0] = real(V);
  I[1] = imag(V);
  I[2] = real(VInv);
  I[3] = imag(VInv);
  I[4] = Sigma; 
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetInvChiDotbB(double X[3], double bB[3],
                    SVTensor *EpsSVT, cdouble Omega,
                    cdouble E[3])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Chi[3][3], InvChi[3][3];
  EpsSVT->Evaluate( Omega, X, Chi);
  Chi[0][0] -= 1.0;
  Chi[1][1] -= 1.0;
  Chi[2][2] -= 1.0;
  Invert3x3Matrix(Chi, InvChi);

  cdouble PreFac = II*ZVAC/Omega;
  E[0]=E[1]=E[2]=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    E[Mu] += PreFac*InvChi[Mu][Nu]*bB[Nu];
}

/***************************************************************/
/* fdim = 20 ***************************************************/
/***************************************************************/
void OverlapIntegrand_PFT(double X[3],
                          double bA[3],     double DivbA,
                          double bB[3],     double DivbB,
                          SVTensor *EpsSVT, cdouble Omega,
                          void *UserData,   double *I)
{
  (void) DivbA; 

  /***************************************************************/
  /* with 'E' defined to be PreFac*InvChi*bB, get E, then get    */
  /* its derivatives by finite-differencing                      */
  /***************************************************************/
  cdouble E[3];
  GetInvChiDotbB(X, bB, EpsSVT, Omega, E);

  cdouble dE[3][3];
  for(int Mu=0; Mu<3; Mu++)
   { 
     double XP[3], DeltaX = (X[Mu]==0.0 ? 1.0e-4 : 1.0e-4*fabs(X[Mu]));
     XP[0] = X[0];
     XP[1] = X[1];
     XP[2] = X[2];
     XP[Mu] += DeltaX;

     cdouble EP[3], EM[3];
     GetInvChiDotbB(XP, bB, EpsSVT, Omega, EP);
     XP[Mu] -= 2.0*DeltaX;
     GetInvChiDotbB(XP, bB, EpsSVT, Omega, EM);

     for(int Nu=0; Nu<3; Nu++)
      dE[Mu][Nu] = (EP[Nu]-EM[Nu])/(2.0*DeltaX);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Chi[3][3], InvChi[3][3];
  EpsSVT->Evaluate( Omega, X, Chi);
  Chi[0][0] -= 1.0;
  Chi[1][1] -= 1.0;
  Chi[2][2] -= 1.0;
  Invert3x3Matrix(Chi, InvChi);
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    dE[Mu][Nu] += (II*ZVAC/Omega) * InvChi[Nu][Mu] * (DivbB/3.0);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double XmX0[3];
  memcpy(XmX0, X, 3*sizeof(double));
  double *XTorque = (double *) UserData;
  if (XTorque) VecPlusEquals(XmX0, -1.0, XTorque);

  cdouble dETheta[3][3];
  for(int Mu=0; Mu<3; Mu++)
   { dETheta[0][Mu] = XTorque[1]*dE[2][Mu] - XTorque[2]*dE[1][Mu];
     dETheta[1][Mu] = XTorque[2]*dE[0][Mu] - XTorque[0]*dE[2][Mu];
     dETheta[2][Mu] = XTorque[0]*dE[1][Mu] - XTorque[1]*dE[0][Mu];
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble *zI = (cdouble *)I;
  memset(zI, 0, 10*sizeof(cdouble));

  double PF1 = 0.5;
  double PF2 = 0.5*TENTHIRDS/real(Omega);
  for(int Mu=0; Mu<3; Mu++)
   { 
     zI[0]   += PF1*bA[Mu]*E[Mu];

     zI[1+0] += PF2*bA[Mu]*dE[0][Mu];
     zI[1+1] += PF2*bA[Mu]*dE[1][Mu];
     zI[1+2] += PF2*bA[Mu]*dE[2][Mu];

     zI[4+0] += PF2*bA[Mu]*dETheta[0][Mu];
     zI[4+1] += PF2*bA[Mu]*dETheta[1][Mu];
     zI[4+2] += PF2*bA[Mu]*dETheta[2][Mu];
   };

  // separate contributions of JxE torque term in zI[7,8,9]
  zI[7+0] = PF2*(bA[1]*E[2] - bA[2]*E[1]);
  zI[7+1] = PF2*(bA[2]*E[0] - bA[0]*E[2]);
  zI[7+2] = PF2*(bA[0]*E[1] - bA[1]*E[0]);
  
  zI[4+0] += zI[7+0];
  zI[4+1] += zI[7+1];
  zI[4+2] += zI[7+2];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOverlapIntegrals(SWGVolume *O, int nt,
                         int nQA, double PreFacA,
                         int iQB, double PreFacB,
                         OverlapIntegrand Integrand,
                         int fdim, void *UserData, 
                         cdouble Omega,
                         int NumPts, double *Integrals)
{
  /***************************************************/
  /***************************************************/
  /***************************************************/
  SWGTet *T     = O->Tets[nt];
  double *QA    = O->Vertices + 3*(T->VI[ nQA ]);
  double *V1    = O->Vertices + 3*(T->VI[ (nQA+1)%4 ]);
  double *V2    = O->Vertices + 3*(T->VI[ (nQA+2)%4 ]);
  double *V3    = O->Vertices + 3*(T->VI[ (nQA+3)%4 ]);
  double *QB    = O->Vertices + 3*iQB;

  double L1[3], L2[3], L3[3];
  VecSub(V1, QA, L1);
  VecSub(V2, QA, L2);
  VecSub(V3, QA, L3);

  SVTensor *EpsSVT= O->SVT;

  /***************************************************/
  /***************************************************/
  /***************************************************/
  double *TetCR = GetTetCR(NumPts);
  memset(Integrals,0,fdim*sizeof(double));
  double *dI = new double[fdim];
  for(int np=0; np<NumPts; np++)
   { 
     double u1=TetCR[4*np + 0];
     double u2=TetCR[4*np + 1];
     double u3=TetCR[4*np + 2];
     double w=(6.0*T->Volume)*TetCR[4*np + 3];

     double x[3], bA[3], bB[3];
     for(int Mu=0; Mu<3; Mu++)
      { bA[Mu] = u1*L1[Mu] + u2*L2[Mu] + u3*L3[Mu];
        x[Mu]  = bA[Mu] + QA[Mu];
        bB[Mu] = x[Mu]  - QB[Mu];
        bA[Mu] *= PreFacA;
        bB[Mu] *= PreFacB;
      };

     Integrand(x, bA, 3.0*PreFacA, bB, 3.0*PreFacB,
               EpsSVT, Omega, UserData, dI);

     for(int nf=0; nf<fdim; nf++)
      Integrals[nf] += w*dI[nf];

   }; 

  delete[] dI;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetOverlapEntries(SWGVolume *O, int nfA,
                      OverlapIntegrand Integrand, int fdim,
                      void *UserData, cdouble Omega,
                      int nfBList[MAXOVERLAP],
                      double Entries[][MAXOVERLAP])
{
  int NNZ=1;
  nfBList[0] = nfA;
  for(int nf=0; nf<fdim; nf++)
   Entries[nf][0]=0.0;

  // loop over the two (positive and negative) tetrahedra of bf #nfA
  SWGFace *FA = O->Faces[nfA];
  double *Integrals=new double[fdim];
  for(int SignA=1; SignA>=-1; SignA-=2)
   {
     int nt         = (SignA==1) ? FA->iPTet : FA->iMTet;
     SWGTet *T      = O->Tets[nt];
     int nQA        = SignA==1 ? FA->PIndex : FA->MIndex;
     double PreFacA = SignA * FA->Area / (3.0*T->Volume);

     // loop over the four faces of tetrahedron #ntA
     for(int ifB=0; ifB<4; ifB++)
      { 
        int nfB = T->FI[ifB];
        if (nfB<0 || nfB >= O->NumInteriorFaces) continue;

        SWGFace *FB    = O->Faces[nfB];
        int SignB      = (nt==FB->iPTet) ? 1 : -1;
        int iQB        = SignB==1 ? FB->iQP : FB->iQM;
        double PreFacB = SignB * FB->Area / (3.0*T->Volume);

        // get contributions of this tet to <b_{nfA}|O|b_{nfB}>
        int NumPts=33; // number of cubature points 
        GetOverlapIntegrals(O, nt, nQA, PreFacA, iQB, PreFacB,
                            Integrand, fdim, UserData, Omega,
                            NumPts, Integrals);

        if (nfB==nfA)
         { for(int nf=0; nf<fdim; nf++) 
            Entries[nf][0]+=Integrals[nf];
         }
        else 
         { 
           nfBList[NNZ] = nfB;
           for(int nf=0; nf<fdim; nf++)
            Entries[nf][NNZ]=Integrals[nf];
           NNZ++;
         };

      };
   };
  delete[] Integrals;

  return NNZ;
 
}

} // namespace buff
