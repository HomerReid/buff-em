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
 * GMatrixElements.cc -- libSWG routines for computing matrix elements
 *                    -- of the dyadic Helmholtz operator and its derivatives
 *                    -- between SWG basis functions
 *
 * homer reid         -- 5/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libSGJC.h"
#include "libscuff.h"
#include "libbuff.h"
#include "TTaylorDuffy.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_OPENMP
#  include <omp.h>
#endif

using namespace scuff;

namespace buff {

#define II cdouble(0.0,1.0)

/***************************************************************/
/* integrand routine and user data structure for computing     */
/* G matrix elements using TetTetInt                           */
/***************************************************************/
typedef struct GMEData
 {
   cdouble k;
   bool NeedDerivatives;
   double XTorque[3];
   int rPower;
 } GMEData;
 
void GMEIntegrand(double *xA, double *bA, double DivbA,
                  double *xB, double *bB, double DivbB,
                  void *UserData, double *I)
{
  GMEData *Data        = (GMEData *)UserData;
  bool NeedDerivatives = Data->NeedDerivatives;
  cdouble k            = Data->k;
  int rPower           = Data->rPower;
  cdouble *zI          = (cdouble *)I;

  double R[3]; 
  R[0] = (xA[0] - xB[0]);
  R[1] = (xA[1] - xB[1]);
  R[2] = (xA[2] - xB[2]);

  double r = sqrt( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );
  if ( r < 1.0e-12 )
   { if (Data->NeedDerivatives)
      memset(zI,0,7*sizeof(double)); 
     else
      zI[0]=0.0;
     return;
   };

  double DotProduct = bA[0]*bB[0] + bA[1]*bB[1] + bA[2]*bB[2];
  cdouble PEFIE = DotProduct - DivbA*DivbB/(k*k);
  cdouble IKR = II*k*r;
  cdouble Phi= exp(IKR)/(4.0*M_PI*r);

  if (rPower!=-10) Phi=pow(r,rPower);

  zI[0] = PEFIE*Phi;

  if (NeedDerivatives)
   { 
     cdouble Psi = (IKR-1.0) * Phi / (r*r);
     if (rPower!=-10) Psi=pow(r,rPower);
     zI[1] = R[0]*PEFIE*Psi;
     zI[2] = R[1]*PEFIE*Psi;
     zI[3] = R[2]*PEFIE*Psi;

     double XmXTorque[3];
     XmXTorque[0] = xA[0] - (Data->XTorque[0]);
     XmXTorque[1] = xA[1] - (Data->XTorque[1]);
     XmXTorque[2] = xA[2] - (Data->XTorque[2]);
     for(int Mu=0; Mu<3; Mu++)
      { int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
        zI[4 + Mu]
         = (XmXTorque[MP1]*R[MP2]-XmXTorque[MP2]*R[MP1])*PEFIE*Psi;
      };
   };

}

/***************************************************************/
/* Use the Taylor-Duffy method to compute the contribution of  */
/* a single tetrahedron-tetrahedron pair to the matrix         */
/* element of G and possible its derivatives.                  */
/* GMETTI stands for 'G matrix element tet-tet integral.'      */
/***************************************************************/
void GetGMETTI_TaylorDuffy(SWGVolume *VA, int OVIA[4], int iQA,
                           SWGVolume *VB, int OVIB[4], int iQB,
                           cdouble k, int ncv, bool *NeedDerivatives,
                           int rPower, cdouble TTI[7])
{
  /*-----------------------------------------------------------*/
  /* initialize Taylor-Duffy argument structure                */
  /*-----------------------------------------------------------*/
  TTDArgStruct MyArgs, *Args=&MyArgs;
  InitTTDArgs(Args);

  if (VA->OTGT) VA->OTGT->Apply(Args->XTorque);
  if (VA->GT) VA->GT->Apply(Args->XTorque);

  int PIndex[14], KIndex[14];
  cdouble KParam[14];
  PIndex[0] = TTD_BDOTBP;
  KIndex[0] = TTD_HELMHOLTZ;
  PIndex[1] = TTD_UNITY;
  KIndex[1] = TTD_HELMHOLTZ;
  int NumPKs = 2;
  for(int Mu=0; Mu<3; Mu++)
   { if (NeedDerivatives && NeedDerivatives[Mu])
      { PIndex[NumPKs]   = TTD_RXBDOTBP + Mu;
        PIndex[NumPKs+1] = TTD_RXUNITY  + Mu;
        KIndex[NumPKs]   = TTD_GRADHELMHOLTZ;
        KIndex[NumPKs+1] = TTD_GRADHELMHOLTZ;
        NumPKs+=2;
      };
     if (NeedDerivatives && NeedDerivatives[3+Mu])
      { PIndex[NumPKs]   = TTD_TXBDOTBP + Mu;
        PIndex[NumPKs+1] = TTD_TXUNITY  + Mu;
        KIndex[NumPKs]   = TTD_GRADHELMHOLTZ;
        KIndex[NumPKs+1] = TTD_GRADHELMHOLTZ;
        NumPKs+=2;
      };
   };
  for(int npk=0; npk<NumPKs; npk++)
   KParam[npk] = k;
 
  if (rPower!=-10)
   { for(int npk=0; npk<NumPKs; npk++)
      { KIndex[npk]=TTD_RP;
        KParam[npk]=rPower;
      };
   };

  Args->NumPKs = NumPKs;
  Args->PIndex = PIndex;
  Args->KIndex = KIndex;
  Args->KParam = KParam;
            
  Args->WhichCase=ncv;
  Args->V1     = VA->Vertices             + 3*OVIA[0];
  Args->V2     = Args->V2P = VA->Vertices + 3*OVIA[1];
  Args->V3     = Args->V3P = VA->Vertices + 3*OVIA[2];
  Args->V4     = Args->V4P = VA->Vertices + 3*OVIA[3];
  if (ncv<4)     Args->V4P = VB->Vertices + 3*OVIB[3];
  if (ncv<3)     Args->V3P = VB->Vertices + 3*OVIB[2];
  if (ncv<2)     Args->V2P = VB->Vertices + 3*OVIB[1];
  Args->Q      = VA->Vertices + 3*iQA;
  Args->QP     = VB->Vertices + 3*iQB;

  cdouble TDI[14], Error[14];
  Args->Result = TDI;
  Args->Error  = Error;

  // calculate taylor-duffy integrals
  TTaylorDuffy(Args);

  // assemble results into output quantities
  cdouble NOK2=9.0/(k*k);
  TTI[0] = TDI[0] - NOK2*TDI[1];

  if (NeedDerivatives)
   { int npk=2;
     for(int Mu=0; Mu<3; Mu++)
      { if (NeedDerivatives[Mu])
        { TTI[1 + Mu] = TDI[npk+0] - NOK2*TDI[npk+1];
          npk+=2;
        };
       if (NeedDerivatives[3+Mu])
        { TTI[4 + Mu] = TDI[npk+0] - NOK2*TDI[npk+1];
          npk+=2;
        };
      };
   };

}

/***************************************************************/
/* If NeedDerivative and dG are nonzero, they must point to    */
/* vectors of length 6. The index Mu into this vector has the  */
/* following significance:                                     */
/*  Mu=0,1,2 -->  d/dx, d/dy, d/dz                             */
/*  Mu=3,4,5 -->  d/dTheta_x, d/dTheta_y, d/dTheta_z           */
/*                                                             */
/* if rPower!=-10, then the usual helmholtz kernel is replaced */
/* with r^rPower.                                              */
/***************************************************************/
cdouble GetGMatrixElement(SWGVolume *VA, int nfA,
                          SWGVolume *VB, int nfB,
                          cdouble Omega, 
                          bool *NeedDerivatives,
                          cdouble *dG, int rPower, bool ForceBF)
{
  cdouble Result[7], Error[7];

  // derivatives vanish identically for diagonal matrix elements
  if (VA==VB && nfA==nfB && NeedDerivatives)
   { if (dG) memset(dG, 0, 6*sizeof(cdouble));
     NeedDerivatives=false; 
   };

  GMEData MyData, *Data=&MyData;
  Data->k=Omega;
  Data->NeedDerivatives = (NeedDerivatives == 0) ? false : true;
  Data->rPower=rPower;
  int fdim = Data->NeedDerivatives ? 7 : 1;


  if (NeedDerivatives) 
   { Data->XTorque[0]=Data->XTorque[1]=Data->XTorque[2]=0.0;
     if (VA->OTGT) VA->OTGT->Apply(Data->XTorque);
     if (VA->GT)   VA->GT->Apply(Data->XTorque);
   };

  double rRel;
  int ncv = CompareBFs(VA, nfA, VB, nfB, &rRel);

  if ( ncv==0 )
   { 
     /***************************************************************/
     /* no common vertices; do a single 6D cubature over the        */
     /* supports of both basis functions                            */
     /***************************************************************/
     int NumPts=4;
     BFBFInt(VA, nfA, VB, nfB, GMEIntegrand, (void *)Data, 2*fdim,
             (double *)Result, (double *)Error, NumPts, 0, 0.0);
    
     if (dG) memcpy(dG, Result+1, 6*sizeof(cdouble));
     return Result[0];
   };

  /***************************************************************/
  /* one or more common vertices; compute the four tet-tet       */
  /* integrals separately                                        */
  /***************************************************************/
  SWGFace *FA = VA->Faces[nfA];
  SWGFace *FB = VB->Faces[nfB];
  double AreaFactor=FA->Area * FB->Area;
  for(int ASign=0; ASign<2; ASign++)
   for(int BSign=0; BSign<2; BSign++)
    { 
      int ntA = (ASign==0) ? FA->iPTet : FA->iMTet;
      int ntB = (BSign==0) ? FB->iPTet : FB->iMTet;
      int OVIA[4], OVIB[4];
      ncv=CompareTets(VA, ntA, VB, ntB, OVIA, OVIB);

      /***************************************************************/
      /* compute the tet--tet integral using brute-force cubature or */
      /* taylor-duffy                                                */
      /***************************************************************/
      cdouble TTI[7]; // tetrahedron-tetrahedron integral
      if (ncv<=1 || ForceBF)
       { 
         int NumPts=16;
         int iQA = (ASign==0) ? FA->PIndex : FA->MIndex;
         int iQB = (BSign==0) ? FB->PIndex : FB->MIndex;
         TetTetInt(VA, ntA, iQA, 1.0, VB, ntB, iQB, 1.0,
                   GMEIntegrand, (void *)Data, 2*fdim,
                   (double *)TTI, (double *)Error, NumPts, 0, 0);
       }
      else
       { 
         int iQA = (ASign==0) ? FA->iQP   : FA->iQM;
         int iQB = (BSign==0) ? FB->iQP   : FB->iQM;
         GetGMETTI_TaylorDuffy(VA, OVIA, iQA, VB, OVIB, iQB,
                               Omega, ncv, NeedDerivatives,
                               rPower, TTI);
         for(int nf=0; nf<fdim; nf++)
          TTI[nf] *= 4.0*AreaFactor;

       }; // if (ncv<=1 ... else )

      double Sign = (ASign==BSign) ? 1.0 : -1.0;
      for(int nf=0; nf<fdim; nf++)
       Result[nf] += Sign*TTI[nf];

    }; // for(int ASign...for(int BSign...)

  if (dG) memcpy(dG, Result+1, 6*sizeof(cdouble));
  return Result[0];

}

} // namespace buff
