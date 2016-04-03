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
 *
 */

#include "libSpherical.h"
#include "buff-neq.h"

#define II cdouble(0.0,1.0)

using namespace buff;

namespace buff {
double GetThetaFactor(double Omega, double T);
               }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSMN(double k, double R, double SMN[2])
{
  double a=k*R, a2=a*a, a3=a2*a, a4=a2*a2, a5=a4*a, a6=a5*a, a7=a6*a, a9=a7*a2;
  double k3=k*k*k;
  if (a<0.1)
   { SMN[0] = (a5/45.0 - a7/315.0 + a9/4725.0) / k3;
     SMN[1] = (2.0*a3/9.0 - 2.0*a5/45.0 + a7/225.0 - 11.0*a9/42525.0) / k3;
   }
  else
   { SMN[0]=fabs(-2 + 2*a2 + 2*cos(2*a) + a*sin(2*a))/(4*a*k*k*k);
     SMN[1]=fabs( (-1 + a4 + cos(2*a) 
                      - a*cos(a)*(2*a*cos(a) + (-4 + a2)*sin(a))
                  )/(2*a3*k3)
                );
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData
 { double Omega;
   int lMax;
   SVTensor *EpsSVT, *TemperatureSVT;
   double ThetaEnvironment;
 } IntegrandData;

void Integrand(double *x, double *b, double Divb,
               void *UserData, double *I)

{ 
  (void )Divb;
  (void )b;

  IntegrandData *Data      = (IntegrandData*)UserData;
  double Omega             = Data->Omega;
  SVTensor *TemperatureSVT = Data->TemperatureSVT;
  int lMax                 = Data->lMax;
  int NAlpha               = (lMax+1)*(lMax+1);
  int NBF                  = 2*(NAlpha-1);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble EpsM1[3][3];
  Data->EpsSVT->Evaluate(Omega, x, EpsM1);
  EpsM1[0][0]-=1.0;
  EpsM1[1][1]-=1.0;
  EpsM1[2][2]-=1.0;

  double T = TemperatureSVT ? TemperatureSVT->EvaluateD(0,x) : 0.0;
  double ThetaEnvironment = Data->ThetaEnvironment;
  double DeltaTheta = GetThetaFactor(Omega, T ) - ThetaEnvironment;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double r, Theta, Phi;
  CoordinateC2S(x, &r, &Theta, &Phi);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (lMax>3) ErrExit("%s:%i:internal error");
#define WSSIZE 20
#define MNSIZE 48
  double Workspace[WSSIZE];
  cdouble MArray[MNSIZE], NArray[MNSIZE], BFArray[2*MNSIZE];
  double k = Data->Omega;
  GetMNlmArray(lMax, k, r,  Theta,  Phi,
               LS_REGULAR, MArray,  NArray, Workspace);

  double SMN[2];
  GetSMN(k, 1.0, SMN);
  double MNorm=1.0/sqrt(SMN[0]);
  double NNorm=1.0/sqrt(SMN[1]);
 
  for(int Alpha=1; Alpha<NAlpha; Alpha++)
   {
     cdouble MC[3], NC[3];
     VectorS2C(Theta, Phi, MArray + 3*Alpha, MC);
     VectorS2C(Theta, Phi, NArray + 3*Alpha, NC);
    
     int nbf = 0 + (Alpha-1);
     BFArray[3*nbf + 0] = MC[0]*MNorm;
     BFArray[3*nbf + 1] = MC[1]*MNorm;
     BFArray[3*nbf + 2] = MC[2]*MNorm;

     nbf = (NAlpha-1) + (Alpha-1);
     BFArray[3*nbf + 0] = NC[0]*NNorm;
     BFArray[3*nbf + 1] = NC[1]*NNorm;
     BFArray[3*nbf + 2] = NC[2]*NNorm;
   };

  cdouble *zI=(cdouble *)I;
  for(int nbfa=0, n=0; nbfa<NBF; nbfa++)
   for(int nbfb=0; nbfb<NBF; nbfb++, n++)
    { 
      cdouble *bA = BFArray + 3*nbfa;
      cdouble *bB = BFArray + 3*nbfb;
      
      zI[n]=0.0;
      for(int Mu=0; Mu<3; Mu++)
       for(int Nu=0; Nu<3; Nu++)
         zI[n] += conj(bA[Mu])*EpsM1[Mu][Nu]*bB[Nu];
      zI[n] *= -k*k;
    };

  for(int n=0; n<NBF*NBF; n++)
   zI[n+NBF*NBF] = -2.0*DeltaTheta*imag(zI[n])/(M_PI*Omega*ZVAC);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetVSWFlux(BNEQData *BNEQD, cdouble Omega, double PF[2])
{
  SWGGeometry *G=BNEQD->G;
  SWGVolume   *O=G->Objects[0];

  /***************************************************************/
  /* the total number of SW indices (\ell,m) with l<=lMax is     */
  /* NAlpha = (\ell+1)^2.                                        */
  /* This includes l==0, which we neglect when computing the     */
  /* V matrix.                                                   */
  /* There are two basis functions for every (\ell,m) pair.      */
  /***************************************************************/
int lMax=1;
  int NAlpha=(lMax+1)*(lMax+1);
  int NBF = 2*(NAlpha-1);
  struct IntegrandData MyData, *Data=&MyData;
  Data->Omega            = real(Omega);
  Data->lMax             = lMax;
  Data->EpsSVT           = O->SVT;
  Data->TemperatureSVT   = BNEQD->TemperatureSVTs[0];
  double TEnvironment    = BNEQD->TEnvironment;
  Data->ThetaEnvironment = GetThetaFactor(real(Omega), TEnvironment);

  static int lMaxSaved=0;
  static HMatrix *V=0, *R=0, *W=0, *D;
  if (lMax!=lMaxSaved)
   { lMaxSaved=lMax;
     V=new HMatrix(NBF, NBF, LHM_COMPLEX);
     R=new HMatrix(NBF, NBF, LHM_COMPLEX);
     W=new HMatrix(NBF, NBF, LHM_COMPLEX);
     D=new HMatrix(NBF, NBF, LHM_COMPLEX);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  V->Zero();
  R->Zero();
  for(int nt=0; nt<O->NumTets; nt++)
   { 
     int fdim=2*NBF*NBF;
#define NBFMAX 30
     cdouble VRElements[2*NBFMAX*NBFMAX], Error[2*NBFMAX*NBFMAX];

     TetInt(O, nt, 0, 1.0, Integrand, (void *)Data,
            2*fdim, (double *)VRElements, (double *)Error, 
            33, 0, 0);
     
     int n=0;
     for(int p=0; p<NBF; p++)
      for(int q=0; q<NBF; q++)
       V->AddEntry(p, q, VRElements[n++]);

     for(int p=0; p<NBF; p++)
      for(int q=0; q<NBF; q++)
       R->AddEntry(p, q, VRElements[n++]);
   };

  /***************************************************************/
  /* fill in the G matrix                                        */
  /***************************************************************/
  double RR=1.0;
  double k=real(Omega), k2=k*k;
  double a=k*RR, a2=a*a, a3=a2*a, a4=a3*a2, a5=a4*a, a6=a5*a, a7=a6*a;

  cdouble GDiag[NBFMAX];
  for(int n=0; n<3; n++)
   GDiag[n]=-1.0/(3.0*k2) + II*(a5/45.0 - a7/315.0)/k2;
  for(int n=0; n<3; n++)
   GDiag[n]=-1.0/(3.0*k2) +(2.0*a3/9.0 - 2.0*a5/45.0)/k2;

  double ImdzGMN = -(a3/9.0 - 2.0*a5/105.0) / sqrt(10.0);

  /***************************************************************/
  /* compute W = (1+V*G)^-1,                                     */
  /* compute D = W*R*W^\dagger                                   */
  /*  D_{ab} = W_{am} * R_{mn} * (W_{bn}^*)                        */
  /***************************************************************/
  for(int p=0; p<NBF; p++)
   for(int q=0; q<NBF; q++)
    W->SetEntry( p, q, (p==q) ? 1.0 : 0.0 
                       +(V->GetEntry(p,q))*GDiag[q]);

  D->Zero();
  for(int p=0; p<NBF; p++)
   for(int q=0; q<NBF; q++)
    for(int m=0; m<NBF; m++)
     for(int n=0; n<NBF; n++)
      D->AddEntry(p, q,        W->GetEntry(p,m) 
                             * R->GetEntry(m,n) 
                         *conj(W->GetEntry(q,n)));

  /***************************************************************/
  /* power = -kz/2 * trace( re(D) * im(G) )                      */
  /* zforce = -z/2c * trace( im(D) * im(dG) )                    */
  /***************************************************************/
  double Power = 0.0, zForce=0.0;
  for(int p=0; p<NBF; p++)
   Power += imag(GDiag[p])*real(D->GetEntry(p,p));
  Power*=0.5*real(Omega)*ZVAC;

  zForce= ZVAC*TENTHIRDS*ImdzGMN*imag( D->GetEntry(0,3) - D->GetEntry(2,6));
  
  PF[0]=Power;
  PF[1]=zForce;

}
