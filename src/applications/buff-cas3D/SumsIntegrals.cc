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
 * SumsIntegrals.cc -- routines for evaluating matsubara sums and
 *                  -- imaginary-frequency integrals
 *
 * homer reid       -- 5/2014
 *
 */

#include "buff-cas3D.h"
#include <stdlib.h>
#include <libSGJC.h>

#define XIMIN  0.001
#define XIMAX 10.000

/***************************************************************/
/* CacheRead: attempt to bypass an entire GetXiIntegrand       */
/* calculation by reading results from the .byXi file.         */
/* Returns 1 if successful (which means the values of the      */
/* energy/force/torque integrand for ALL transformations at    */
/* this value of Xi were successfully read from the file) or 0 */
/* on failure).                                                */
/***************************************************************/
int CacheRead(const char *ByXiFileName, BC3Data *BC3D, 
              double Xi, double *EFT)
{ 
  if (BC3D->UseExistingData==false)
   return 0;

  FILE *f;
  double Q[4];
  char Line[1000], fTag[1000];
  double fXi;
  int nt, ntnq, nRead, FoundFirst;

  /*----------------------------------------------------------*/
  /* 0. try to open the cache file. --------------------------*/
  /*----------------------------------------------------------*/
  if ( !(f=fopen(ByXiFileName,"r")) )
   return 0;

  /*----------------------------------------------------------*/
  /*----------------------------------------------------------*/
  /*----------------------------------------------------------*/
  for(;;)
   {
     /*----------------------------------------------------------*/
     /* 1. skip down through the cache file until we find a line */
     /*    whose Xi and Tag values equal Xi and the first        */
     /*    tag in the workspace structure.                       */
     /*----------------------------------------------------------*/
     FoundFirst=0;
     while( !FoundFirst && fgets(Line,1000,f) )
      { sscanf(Line,"%s %le",fTag,&fXi);
        if ( fabs(fXi-Xi) < 1.0e-8*Xi && !strcmp(fTag,BC3D->GTCList[0]->Tag) )
         FoundFirst=1;
      };
     if ( !FoundFirst ) 
      { fclose(f); 
        return 0;
      };
   
     Log(" found (Tag,Xi)=(%s,%e) in cache file...",fTag,Xi);
   
     /*----------------------------------------------------------*/
     /* 2. verify that the line we just read from the cache file */
     /*    contains data for all the quantities we need          */
     /*----------------------------------------------------------*/ 
     nRead=sscanf(Line,"%s %le %le %le %le %le",fTag,&fXi,Q,Q+1,Q+2,Q+3); 
     if ( nRead != BC3D->NumQuantities+2 )
      { Log(" ...but number of quantities is wrong (skipping)");
        continue;
      };
     memcpy(EFT,Q,BC3D->NumQuantities*sizeof(double));
     ntnq=BC3D->NumQuantities;
   
     /*----------------------------------------------------------*/
     /* 3. ok, since that worked, now keep going ----------------*/
     /*----------------------------------------------------------*/
     for(nt=1; nt<BC3D->NumTransformations; nt++)
      { 
        /* check for premature end of file */
        if ( !fgets(Line,1000,f) )
         { Log(" ...but data for some transforms were missing (skipping)");
           break;
         };
   
        nRead=sscanf(Line,"%s %le %le %le %le %le",fTag,&fXi,Q,Q+1,Q+2,Q+3);
   
        /* check for incorrect number of quantities */
        if ( nRead != BC3D->NumQuantities+2 )
         { Log(" ...but number of quantities is wrong (skipping)");
           break;
         };
   
        /* check for tag and/or Xi mismatch */
        if ( fabs(fXi-Xi)>1.0e-8*Xi || strcmp(fTag,BC3D->GTCList[nt]->Tag) )
         { Log(" ...but tag #%i did not match (%s != %s) (skipping)",
               nt,fTag,BC3D->GTCList[nt]->Tag);
           break;
         };
   
        memcpy(EFT+ntnq,Q,BC3D->NumQuantities*sizeof(double));
        ntnq+=BC3D->NumQuantities;
      };
   
     if (ntnq==BC3D->NTNQ)
      { Log(" ...and successfully read data for all quantities at all transforms");
        fclose(f);
        return 1;
      };

   };

}

/***************************************************************/
/* wrapper with correct prototype for pcubature                */
/***************************************************************/
int GetXiIntegrand_Adaptive(unsigned ndim, const double *x, 
                            void *params, unsigned fdim, double *fval)
{
  (void) ndim; // unused
  (void) fdim; // unused

  if (x[0]==1.0)
   { memset(fval, 0, fdim*sizeof(double));
     return 0;
   };

  double Xi = XIMIN + x[0]/(1.0-x[0]);
  double Jacobian = 1.0/( (1.0-x[0])*(1.0-x[0]) );
  BC3Data *BC3D = (BC3Data *)params;
  double *EFT = fval;

  GetXiIntegrand(BC3D, Xi, EFT);

  for(int ntnq=0; ntnq<BC3D->NTNQ; ntnq++)
   EFT[ntnq]*=Jacobian;

  return 0;

}

/***************************************************************/
/* Evaluate the Matsubara sum to get total Casimir quantities  */
/* at temperature Temperature degrees Kelvin.                  */
/*                                                             */
/* how the temperature conversion works:                       */
/*  a. temperature in eV = kT = 8.6173e-5 * (T in kelvin)      */
/*  b. temperature in our internal energy units                */
/*     = (kT in eV) / (0.1973 eV)                              */
/*     = (8.6173e-5 / 0.1973 ) * (T in Kelvin)                 */
/***************************************************************/
#define BOLTZMANNK 4.36763e-4
void GetMatsubaraSum(BC3Data *BC3D, double Temperature, double *EFT, double *Error)
{ 
  int n, ntnq, NTNQ=BC3D->NTNQ;
  double Xi, Weight;

  double *dEFT = new double [NTNQ]; 
  double *LastEFT = new double[NTNQ]; 
  double RelDelta;
  int AllConverged=0;
  int *ConvergedIters = new int [NTNQ];

  double kT = BOLTZMANNK * Temperature;

  memset(EFT,0,NTNQ*sizeof(double));
  memset(BC3D->Converged,0,NTNQ*sizeof(int));
  memset(ConvergedIters,0,NTNQ*sizeof(int));

  Log("Beginning Matsubara sum at T=%g kelvin...",Temperature);

  for(n=0; n<BC3D->MaxXiPoints; n++)
   { 
     /***************************************************************/
     /* compute the next matsubara frequency ************************/
     /***************************************************************/
     if (n==0)
      { Weight=0.5;
        // NOTE: we assume that the integrand is constant for Xi < XIMIN
        Xi=XIMIN;
      }
     else
      { Weight=1.0;
        Xi=2.0*M_PI*kT*((double)n);
      };

     /***************************************************************/
     /* evaluate the frequency integrand at this matsubara frequency*/
     /***************************************************************/
     GetXiIntegrand(BC3D, Xi, dEFT);

     /***************************************************************/
     /* accumulate contributions to the sum.                        */
     /*                                                             */
     /* how it works: the matsubara sum is                          */
     /*  2\pi kT *  \sum_n^\prime F(\xi_n)                          */
     /* where \xi_n is the nth matsubara frequency and F(\xi) is    */
     /* the casimir integrand (and the primed sum means that the    */
     /* n==0 term is weighted with a factor of 1/2).                */
     /*                                                             */
     /* however, my GetXiIntegrand() routine returns the quantity   */
     /* FI = F(\xi_n) / (2\pi). (this is so that the integral of FI */
     /* over all \xi returns the correct casimir quantity with no   */
     /* additional multiplicative prefactors.)                      */
     /*                                                             */
     /* thus the matsubara sum is                                   */
     /*  4\pi^2 kT *  \sum_n^\prime FI(\xi_n)                       */
     /*                                                             */
     /* where FI is what is returned by GetXiIntegrand.             */
     /***************************************************************/
     memcpy(LastEFT,EFT,NTNQ*sizeof(double));
     for(ntnq=0; ntnq<NTNQ; ntnq++)
      EFT[ntnq] += Weight * 4.0*M_PI*M_PI* kT * dEFT[ntnq];

     /*********************************************************************/
     /* convergence analysis.                                             */
     /* how it works: if the relative change in output quantity #ntnq is  */
     /* less than EPSREL, we increment ConvergedIters[ntnq]; otherwise we */
     /* set ConvergedIters[ntnq] to 0.                                    */
     /* when ConvergedIters[ntnq] hits 2, we mark quantity #ntnq as       */
     /* having converged.                                                 */
     /* when all quantities have converged, we are done.                  */
     /*********************************************************************/
     for(AllConverged=1, ntnq=0; ntnq<NTNQ; ntnq++)
      { 
        if ( BC3D->Converged[ntnq] ) 
         continue;

        Error[ntnq] = fabs(EFT[ntnq] - LastEFT[ntnq]);
        RelDelta = Error[ntnq] / fabs(EFT[ntnq]);
        if ( RelDelta < 1.0e-6 )
         ConvergedIters[ntnq]++;
        else
         ConvergedIters[ntnq]=0;

        if ( ConvergedIters[ntnq]>=2 )
         BC3D->Converged[ntnq]=1;
        else  
         AllConverged=0;
      }; 

     if (AllConverged==1)
      break;

   }; /* for (n=0 ... */

  delete[] dEFT;   
  delete[] LastEFT;
  delete[] ConvergedIters;
  
  if (AllConverged==0)
   { 
     fprintf(stderr,"\n*\n* WARNING: Matsubara sum unconverged after %i frequency samples.\n*\n",n);
     Log("Matsubara sum UNCONVERGED at n=%i samples",n); 
   } 
  else 
   Log("Matsubara sum converged after summing n=%i frequency points.",n);
    
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetXiIntegral_Fixed(BC3Data *BC3D, int NumIntervals, double *I, double *E)
{ 
  int fdim = BC3D->NTNQ;
  double *fLeft  = new double[fdim];
  double *fMid   = new double[fdim];
  double *fRight = new double[fdim];
  double Xi;

  /*--------------------------------------------------------------*/
  /*- evaluate integrand at leftmost frequency point and estimate */
  /*- the integral from 0 to XIMIN by assuming that the integrand */
  /*- is constant in that range                                   */
  /*--------------------------------------------------------------*/
  Xi=XIMIN;
  GetXiIntegrand(BC3D, Xi, fLeft);
  for(int nf=0; nf<fdim; nf++)
   I[nf] = fLeft[nf] * XIMIN;
  memset(E,0,BC3D->NTNQ*sizeof(double));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Delta = (XIMAX - XIMIN) / NumIntervals;
  double *ISimp  = new double[fdim];
  double *ITrap  = new double[fdim];
  for(int nIntervals=0; nIntervals<NumIntervals; nIntervals++)
   { 
     // evaluate integrand at midpoint of interval 
     Xi += 0.5*Delta;
     GetXiIntegrand(BC3D, Xi, fMid);

     // evaluate integrand at right end of interval 
     Xi += 0.5*Delta;
     GetXiIntegrand(BC3D, Xi, fRight);

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

}

/***************************************************************/
/* Integrate over the positive imaginary frequency axis to get */
/* the total Casimir quantities at zero temperature.           */
/***************************************************************/
void GetXiIntegral_Adaptive(BC3Data *BC3D, double *EFT, double *Error)
{
  double Lower[1] = {0.0}; 
  double Upper[1] = {1.0};

  pcubature(BC3D->NTNQ, GetXiIntegrand_Adaptive, (void *)BC3D, 
            1, Lower, Upper, BC3D->MaxXiPoints, BC3D->AbsTol, 
            BC3D->RelTol, ERROR_INDIVIDUAL, EFT, Error);
   
}
