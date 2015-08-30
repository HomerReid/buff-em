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
 * FrequencyIntegral.cc -- buff-neq module for numerical quadrature 
 *                      -- over frequencies
 *
 * homer reid           -- 5/2012
 *
 */

#include <stdlib.h>
#include <libSGJC.h>
#include "buff-neq.h"

// how this works: 
//  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)
//  b. temperature in our internal energy units
//     = (kT in eV) / (0.1973 eV)
//     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)
#define BOLTZMANNK 4.36763e-4

// for example: suppose in the real world we have 
// omega=3e14 rad/sec, T = 300 kelvin. then
// \hbar \omega / (kT) 
//   = (6.6e-16 ev s )(3e14 s^{-1}) / (0.026 ev) 
//   = 7.6
// whereas in this code we would have Omega==1, T=300, and hence
// Omega/(BOLTZMANNK*T) = (1/(4.36763e-4*300)) = 7.6. 

/***************************************************************/
/***************************************************************/
/***************************************************************/
double Theta(double Omega, double T)
{ 
  if (T==0.0)
   return 0.0;
  return Omega / ( exp( Omega/(BOLTZMANNK*T) ) - 1.0 );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PutInThetaFactors(BNEQData *BNEQD, double Omega,
                       double *TObjects, double TEnvironment,
                       double *FluxVector)
{
  /*--------------------------------------------------------------*/
  /*- quantities arising from sources inside object nss are       */
  /*- weighted by a factor [Theta(T) - Theta(TEnv)]               */
  /*- note: nss = 'num surface, source'                           */
  /*-       nsd = 'num surface, destination'                      */
  /*--------------------------------------------------------------*/
  int NO = BNEQD->G->NumObjects;
  int NT = BNEQD->NumTransformations;
  int NQ = BNEQD->NQ;
  for(int nss=0; nss<NO; nss++)
   { double DeltaTheta = Theta(Omega, TObjects[nss]) - Theta(Omega, TEnvironment);
     for(int nt=0; nt<NT; nt++)
      for(int nsd=0; nsd<NO; nsd++)
       for(int nq=0; nq<NQ; nq++)
        FluxVector[ GetIndex(BNEQD, nt, nss, nsd, nq) ]*= DeltaTheta/M_PI;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=vfopen("%s.integrand","a",BNEQD->FileBase);
  for(int nt=0; nt<NT; nt++)
   for(int nss=0; nss<NO; nss++)
    for(int nsd=0; nsd<NO; nsd++)
     { 
       fprintf(f,"%s %e %i%i ",BNEQD->GTCList[nt]->Tag,Omega,nss+1,nsd+1);
       for(int nq=0; nq<NQ; nq++)
        fprintf(f,"%.8e ",FluxVector[ GetIndex(BNEQD, nt, nss, nsd, nq) ] );
       fprintf(f,"\n");
     };
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteDataToOutputFile(BNEQData *BNEQD, double *I, double *E)
{
  time_t MyTime;
  struct tm *MyTm;
  char TimeString[30];
  
  MyTime=time(0);
  MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  FILE *f=vfopen("%s.out","a",BNEQD->FileBase);
  fprintf(f,"\n");
  fprintf(f,"# buff-neq run on %s (%s)\n",GetHostName(),TimeString);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 transform tag\n");
  fprintf(f,"# 2 (sourceObject, destObject) \n");
  int nq=3;
  if (BNEQD->QuantityFlags & QFLAG_PABS)
   { fprintf(f,"# (%i,%i) power (value,error)\n",nq,nq+1); nq+=2; }
  if (BNEQD->QuantityFlags & QFLAG_XFORCE) 
   { fprintf(f,"# (%i,%i) x-force (value,error)\n",nq, nq+1); nq+=2; }
  if (BNEQD->QuantityFlags & QFLAG_YFORCE) 
   { fprintf(f,"# (%i,%i) y-force (value,error)\n",nq, nq+1); nq+=2; }
  if (BNEQD->QuantityFlags & QFLAG_ZFORCE) 
   { fprintf(f,"# (%i,%i) z-force (value,error)\n",nq, nq+1); nq+=2; }
  if (BNEQD->QuantityFlags & QFLAG_XTORQUE) 
   { fprintf(f,"# (%i,%i) x-torque (value,error)\n",nq, nq+1); nq+=2; }
  if (BNEQD->QuantityFlags & QFLAG_YTORQUE) 
   { fprintf(f,"# (%i,%i) y-torque (value,error)\n",nq, nq+1); nq+=2; }
  if (BNEQD->QuantityFlags & QFLAG_ZTORQUE) 
   { fprintf(f,"# (%i,%i) z-torque (value,error)\n",nq, nq+1); nq+=2; }

  int NO = BNEQD->G->NumObjects;
  int NT = BNEQD->NumTransformations;
  int NQ = BNEQD->NQ;
  double TotalQuantity[MAXQUANTITIES], TotalError[MAXQUANTITIES];
  for(int nt=0; nt<NT; nt++)
   for(int nsd=0; nsd<NO; nsd++)
    { 
      memset(TotalQuantity,0,NQ*sizeof(double));
      memset(TotalError,   0,NQ*sizeof(double));
      for(int nss=0; nss<NO; nss++)
       { fprintf(f,"%s %i%i ",BNEQD->GTCList[nt]->Tag,nss+1,nsd+1);
         for(nq=0; nq<NQ; nq++)
          { int i = GetIndex(BNEQD, nt, nss, nsd, nq);
            fprintf(f,"%+16.8e %+16.8e ", I[i], E[i] );
            TotalQuantity[nq] += I[i];
            TotalError[nq] += E[i];
          };
         fprintf(f,"\n");
       };

      fprintf(f,"%s 0%i ",BNEQD->GTCList[nt]->Tag,nsd+1);
      for(nq=0; nq<NQ; nq++)
       fprintf(f,"%e %e ",TotalQuantity[nq],TotalError[nq]);
      fprintf(f,"\n");

    };
  fclose(f);   
}

#if 0
/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct FIData 
 {
   BNEQData *BNEQD;
   double OmegaMin;
   int Infinite;
   double *TObjects;
   double TEnvironment;

 } FIData;

int SGJCIntegrand(unsigned ndim, const double *x, void *params,
                   unsigned fdim, double *fval)
{
  (void) ndim; // unused
  (void) fdim;

  FIData *FID         = (FIData *)params;

  BNEQData *BNEQD     = FID->BNEQD; 
  double OmegaMin     = FID->OmegaMin;
  int Infinite        = FID->Infinite;
  double *TObjects    = FID->TObjects;
  double TEnvironment = FID->TEnvironment;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Omega, Jacobian; 
  if(Infinite)
   { Omega    = OmegaMin + x[0] / (1.0-x[0]);
     Jacobian = 1.0 / ( (1.0-x[0]) * (1.0-x[0]) );
   }
  else
   { Omega    = x[0];
     Jacobian = 1.0;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GetFlux(BNEQD, Omega, fval);
  for(unsigned int nf=0; nf<fdim; nf++)
   fval[nf]*=Jacobian;
  PutInThetaFactors(BNEQD, Omega, TObjects, TEnvironment, fval);

  return 0.0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateFrequencyIntegral(BNEQData *BNEQD,
                               double OmegaMin, double OmegaMax,
                               double *TObjects, double TEnvironment,
                               double AbsTol, double RelTol,
                               double *I, double *E)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FIData MyFIData, *FID=&MyFIData;
  FID->BNEQD = BNEQD;
  FID->OmegaMin=OmegaMin;

  if (OmegaMax==-1.0)
   { FID->Infinite=1; 
     OmegaMin=0.0;
     OmegaMax=1.0;
   }
  else
   FID->Infinite=0;

  FID->TObjects=TObjects;
  FID->TEnvironment=TEnvironment;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = BNEQD -> G;
  int NO = G->NumObjects;
  int NT = BNEQD->NumTransformations;
  int NQ = BNEQD->NQ;
  int fdim = NT*NO*NO*NQ;
  pcubature_log(fdim, SGJCIntegrand, (void *)FID, 1, &OmegaMin, &OmegaMax,
                1000, AbsTol, RelTol, ERROR_INDIVIDUAL, I, E, "buff-neq.SGJClog");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  WriteDataToOutputFile(BNEQD, I, E);

}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EvaluateFrequencyIntegral2(BNEQData *BNEQD, double OmegaMin, double OmegaMax,
                                double *TObjects, double TEnvironment, int NumIntervals,
                                double *I, double *E)
{ 
  int NO = BNEQD->G->NumObjects;
  int NT = BNEQD->NumTransformations;
  int NQ = BNEQD->NQ;
  int fdim = NT*NO*NO*NQ;
  double *fLeft  = new double[fdim];
  double *fMid   = new double[fdim];
  double *fRight = new double[fdim];
  double Omega;

  if ( OmegaMax == -1.0 )
   { 
     // if the user didn't specify an upper frequency bound, we choose 
     // OmegaMax to be the frequency at which Theta(OmegaMax, TMax) has 
     // decayed to 10^{-10} of the value of Theta(0,TMax), where TMax 
     // is the largest temperature of any object. 
     // note: the function x/(exp(x)-1) falls below 10^{-10} at x=26.3.

     double TMax=TEnvironment;
     for (int ns=0; ns<NO; ns++)
      if ( TObjects[ns]>=0.0 && TObjects[ns]>TMax ) 
       TMax=TObjects[ns];

     if (TMax==0.0) // all temperatures are zero; no point in calculating
      { memset(I,0,fdim*sizeof(double));
        memset(E,0,fdim*sizeof(double));
        return;
      };  

     OmegaMax = 26.3*BOLTZMANNK*TMax;
     Log("Integrating to a maximum frequency of k=%e um^{-1} (w=%e rad/sec)\n",OmegaMax,OmegaMax*3.0e14);

   };

  /*--------------------------------------------------------------*/
  /*- evaluate integrand at leftmost point.  ---------------------*/
  /*- If the leftmost point is less than OMEGAMIN (the smallest  -*/
  /*- frequency at which we can do reliable calculations) then   -*/
  /*- we estimate the integral from OmegaMin to OMEGAMIN by      -*/
  /*- a rectangular rule with integrand values computed at       -*/
  /*- OMEGAMIN.                                                  -*/
  /*--------------------------------------------------------------*/
  #define OMEGAMIN 0.01
  if (OmegaMin < OMEGAMIN)
   { 
     Omega=OMEGAMIN;
     GetFlux(BNEQD, Omega, fLeft);
     PutInThetaFactors(BNEQD, Omega, TObjects, TEnvironment, fLeft);
     for(int nf=0; nf<fdim; nf++)
      { I[nf] = fLeft[nf] * (OMEGAMIN-OmegaMin);
        E[nf] = 0.0;
      };
     OmegaMin=OMEGAMIN;
   }
  else
   { 
     Omega=OmegaMin;
     GetFlux(BNEQD, Omega, fLeft);
     PutInThetaFactors(BNEQD, Omega, TObjects, TEnvironment, fLeft);
     memset(I, 0, fdim*sizeof(double));
     memset(E, 0, fdim*sizeof(double));
   };


  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Delta = (OmegaMax - OmegaMin) / NumIntervals;
  double *ISimp  = new double[fdim];
  double *ITrap  = new double[fdim];
  for(int nIntervals=0; nIntervals<NumIntervals; nIntervals++)
   { 
     // evaluate integrand at midpoint of interval 
     Omega += 0.5*Delta;
     GetFlux(BNEQD, Omega, fMid);
     PutInThetaFactors(BNEQD, Omega, TObjects, TEnvironment, fMid);

     // evaluate integrand at right end of interval 
     Omega += 0.5*Delta;
     GetFlux(BNEQD, Omega, fRight);
     PutInThetaFactors(BNEQD, Omega, TObjects, TEnvironment, fRight);

     // compute the simpson's rule and trapezoidal rule
     // estimates of the integral over this interval  
     // and take their difference as the error
     for(int nf=0; nf<fdim; nf++)
      { ISimp[nf] = (fLeft[nf] + 4.0*fMid[nf] + fRight[nf])*Delta/6.0;
        ITrap[nf] = (fLeft[nf] + 2.0*fMid[nf] + fRight[nf])*Delta/4.0;
        I[nf] += ISimp[nf];
        E[nf] += fabs(ISimp[nf] - ITrap[nf]);
      };

     // prepare for next iteration
     memcpy(fLeft, fRight, fdim*sizeof(double));

   };
  delete[] fLeft;
  delete[] fMid;
  delete[] fRight;
  delete[] ISimp;
  delete[] ITrap;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  WriteDataToOutputFile(BNEQD, I, E);

}
