int RegionOnly=-1;
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
 * TetTetIntegrals_TD.cc -- Taylor-Duffy computation of 
 *                          integrals between pairs of tetrahedra
 *
 * homer reid   -- 5/2014
 */
#include <stdio.h>

#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libbuff.h>
#include <libSGJC.h>
#include <libTriInt.h>

using namespace scuff;

namespace buff{

#define NUMREGIONS 18

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct TDIntegrandData
 { 
   double *V0, L1[3], L2[3], L3[3];
   double *V0P, L1P[3], L2P[3], L3P[3];
   double *Q, *QP, PreFac, PreFacP;

   UserTTIntegrand UserFunction;
   void *UserData;
   double *fBuffer;

 } TDIntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSubregionData(double y1, double y2, double y3,
                      double u[NUMREGIONS][3],
                      double ABCDEF[NUMREGIONS][6],
                      double J[NUMREGIONS])
{
  double y1y2=y1*y2;
  double y1y3=y1*y3;
  double y1y2y3=y1*y2*y3;
  double y12=y1*y1;
  double y12y2=y1*y1*y2;
  double u1, u2, u3;

  // d=0: A, C1, G1 
  u[0][0]=u1=y1y2y3;
  u[0][1]=u2=y1y2;
  u[0][2]=u3=y1;
  ABCDEF[0][0]=-u1 + u3;
  ABCDEF[0][1]=1.0 - u1;
  ABCDEF[0][2]=-u2 + u3;
  ABCDEF[0][3]=u1 - u2;
  ABCDEF[0][4]=0.0;
  ABCDEF[0][5]=u2 - u3;
  J[0]=y12y2;

  // d=1: A, C1, G2
  u[1][0]=u1=y1y2;
  u[1][1]=u2=y1;
  u[1][2]=u3=y1y3;
  ABCDEF[1][0]=-u1 + u2;
  ABCDEF[1][1]=1.0 - u1;
  ABCDEF[1][2]=0.0;
  ABCDEF[1][3]=u1 - u2;
  ABCDEF[1][4]=0.0;
  ABCDEF[1][5]=0.0;
  J[1]=y12;

  // d=2: A, C1, H
  u[2][0]=u1=y1y2y3;
  u[2][1]=u2=y1y2;
  u[2][2]=u3=y1-1.0;
  ABCDEF[2][0]=-u1+u2-u3;
  ABCDEF[2][1]=1.0 - u1;
  ABCDEF[2][2]=-u3;
  ABCDEF[2][3]=u1 - u2;
  ABCDEF[2][4]=-u3;
  ABCDEF[2][5]=0.0;
  J[2]=y12y2;

  // d=3: A, C2, G1
  u[3][0]=u1=y1y2;
  u[3][1]=u2=y1y2y3;
  u[3][2]=u3=1.0+y1y2y3-y1;
  ABCDEF[3][0]=-u2+u3;
  ABCDEF[3][1]=1.0 - u1;
  ABCDEF[3][2]=-u2+u3;
  ABCDEF[3][3]=0.0;
  ABCDEF[3][4]=0.0;
  ABCDEF[3][5]=u2-u3;
  J[3]=y12y2;

  // d=4: A, C2, G2
  u[4][0]=u1=y1;
  u[4][1]=u2=y1y2;
  u[4][2]=u3=y1y2y3;
  ABCDEF[4][0]=0.0;
  ABCDEF[4][1]=1.0 - u1;
  ABCDEF[4][2]=0.0;
  ABCDEF[4][3]=0.0;
  ABCDEF[4][4]=0.0;
  ABCDEF[4][5]=0.0;
  J[4]=y12y2;

  // d=5: A, C2, H
  u[5][0]=u1=y1y2;
  u[5][1]=u2=y1y2y3;
  u[5][2]=u3=y1-1.0;
  ABCDEF[5][0]=-u3;
  ABCDEF[5][1]=1.0 - u1;
  ABCDEF[5][2]=-u3;
  ABCDEF[5][3]=0.0;
  ABCDEF[5][4]=-u3;
  ABCDEF[5][5]=0.0;
  J[5]=y12y2;

  // d=6: A, D, I
  u[6][0]=u1=y1y2y3;
  u[6][1]=u2=y1-1.0;
  u[6][2]=u3=y1y2*(1.0-y3);
  ABCDEF[6][0]=-u2+u3;
  ABCDEF[6][1]=1.0 - u1;
  ABCDEF[6][2]=-u2+u3;
  ABCDEF[6][3]=0.0;
  ABCDEF[6][4]=0.0;
  ABCDEF[6][5]=u2-u3;
  J[6]=y12y2;

  // d=7: A, D, J1
  u[7][0]=u1=y1y2y3;
  u[7][1]=u2=y1y2-1.0;
  u[7][2]=u3=y1-1.0;
  ABCDEF[7][0]=-u2;
  ABCDEF[7][1]=1.0 - u1;
  ABCDEF[7][2]=-u2;
  ABCDEF[7][3]=0.0;
  ABCDEF[7][4]=-u3;
  ABCDEF[7][5]=u2-u3;
  J[7]=y12y2;

  // d=8: A, D, J2
  u[8][0]=u1=y1y2y3;
  u[8][1]=u2=y1-1.0;
  u[8][2]=u3=y1y2-1.0;
  ABCDEF[8][0]=-u3;
  ABCDEF[8][1]=1.0-u1;
  ABCDEF[8][2]=-u3;
  ABCDEF[8][3]=0.0;
  ABCDEF[8][4]=-u3;
  ABCDEF[8][5]=0.0;
  J[8]=y12y2;

  // d=9: B, E, G1
  u[9][0]=u1=y1-1.0;
  u[9][1]=u2=y1y2y3;
  u[9][2]=u3=y1y2;
  ABCDEF[9][0]=-u1+u3;
  ABCDEF[9][1]=1.0;
  ABCDEF[9][2]=-u2+u3;
  ABCDEF[9][3]=u1-u2;
  ABCDEF[9][4]=0.0;
  ABCDEF[9][5]=u2-u3;
  J[9]=y12y2;

  // d=10: B, E, G2
  u[10][0]=u1=y1-1.0;
  u[10][1]=u2=y1y2;
  u[10][2]=u3=y1y2y3;
  ABCDEF[10][0]=-u1+u2;
  ABCDEF[10][1]=1.0;
  ABCDEF[10][2]=0.0;
  ABCDEF[10][3]=u1-u2;
  ABCDEF[10][4]=0.0;
  ABCDEF[10][5]=0.0;
  J[10]=y12y2;

  // d=11: B, E, H
  u[11][0]=u1=y1-1.0;
  u[11][1]=u2=y1*(1.0-y2);
  u[11][2]=u3=-y1y2y3;
  ABCDEF[11][0]=-u1+u2-u3;
  ABCDEF[11][1]=1.0;
  ABCDEF[11][2]=-u3;
  ABCDEF[11][3]=u1-u2;
  ABCDEF[11][4]=-u3;
  ABCDEF[11][5]=0.0;
  J[11]=y12y2;

  // d=12: B, F1, I
  u[12][0]=u1=y1y2-1.0;
  u[12][1]=u2=y1-1.0;
  u[12][2]=u3=y1y2y3;
  ABCDEF[12][0]=-u1+u3;
  ABCDEF[12][1]=1.0;
  ABCDEF[12][2]=-u2+u3;
  ABCDEF[12][3]=u1-u2;
  ABCDEF[12][4]=0.0;
  ABCDEF[12][5]=u2-u3;
  J[12]=y12y2;

  // d=13: B, F1, J1
  u[13][0]=u1=-y1;
  u[13][1]=u2=-y1y2;
  u[13][2]=u3=-y1y2y3;
  ABCDEF[13][0]=-u1;
  ABCDEF[13][1]=1.0;
  ABCDEF[13][2]=-u2;
  ABCDEF[13][3]=u1-u2;
  ABCDEF[13][4]=-u3;
  ABCDEF[13][5]=u2-u3;
  J[13]=y12y2;

  // d=14: B, F1, J2
  u[14][0]=u1=y1y2-1.0; 
  u[14][1]=u2=y1-1.0;
  u[14][2]=u3=y1-1.0-y1y2y3;
  ABCDEF[14][0]=-u1+u2-u3;
  ABCDEF[14][1]=1.0;
  ABCDEF[14][2]=-u3;
  ABCDEF[14][3]=u1-u2;
  ABCDEF[14][4]=-u3;
  ABCDEF[14][5]=0.0;
  J[14]=y12y2;

  // d=15: B, F2, I
  u[15][0]=u1=y1-1.0;
  u[15][1]=u2=y1y2-1.0;
  u[15][2]=u3=y1y2y3;
  ABCDEF[15][0]=-u2+u3;
  ABCDEF[15][1]=1.0;
  ABCDEF[15][2]=-u2+u3;
  ABCDEF[15][3]=0.0;
  ABCDEF[15][4]=0.0;
  ABCDEF[15][5]=u2-u3;
  J[15]=y12y2;

  // d=16: B, F2, J1
  u[16][0]=u1=-y1y2;
  u[16][1]=u2=-y1;
  u[16][2]=u3=-y1y3;
  ABCDEF[16][0]=-u2;
  ABCDEF[16][1]=1.0;
  ABCDEF[16][2]=-u2;
  ABCDEF[16][3]=0.0;
  ABCDEF[16][4]=-u3;
  ABCDEF[16][5]=u2-u3;
  J[16]=y12;

  // d=17: B, F2, J2
  u[17][0]=u1=y1-1.0;
  u[17][1]=u2=y1y2-1.0;
  u[17][2]=u3=y1y2y3-1.0;
  ABCDEF[17][0]=-u3;
  ABCDEF[17][1]=1.0;
  ABCDEF[17][2]=-u3;
  ABCDEF[17][3]=0.0;
  ABCDEF[17][4]=-u3;
  ABCDEF[17][5]=0.0;
  J[17]=y12y2;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetxIntegrand(const double x[3], const double u[3], const double ABCDEF[6],
                   TDIntegrandData *TDIData, double *xJacobian, double *f)
{
  double A  = ABCDEF[0];
  double B  = ABCDEF[1];
  double C  = ABCDEF[2];
  double E  = ABCDEF[4];
  double BmA = B-A;

  double Xi1 = (BmA)*x[0]           + A;
  double Xi2 = (BmA)*x[0]*x[1]      + C;
  double Xi3 = (BmA)*x[0]*x[1]*x[2] + E;
  *xJacobian = BmA*BmA*BmA*x[0]*x[0]*x[1];

  double Eta1=Xi1+u[0];
  double Eta2=Xi2+u[1];
  double Eta3=Xi3+u[2];

  double *V0     = TDIData->V0;
  double *L1     = TDIData->L1;
  double *L2     = TDIData->L2;
  double *L3     = TDIData->L3;
  double *V0P    = TDIData->V0P;
  double *L1P    = TDIData->L1P;
  double *L2P    = TDIData->L2P;
  double *L3P    = TDIData->L3P;
  double *Q      = TDIData->Q;
  double *QP     = TDIData->QP;
  double PreFac  = TDIData->PreFac;
  double PreFacP = TDIData->PreFacP;

  double X[3], b[3], XP[3], bP[3];
  for(int Mu=0; Mu<3; Mu++)
   { X[Mu] = V0[Mu] + Xi1*L1[Mu] + Xi2*L2[Mu] + Xi3*L3[Mu];
     b[Mu] = PreFac * ( X[Mu] - Q[Mu] ); 
     XP[Mu] = V0P[Mu] + Eta1*L1P[Mu] + Eta2*L2P[Mu] + Eta3*L3P[Mu];
     bP[Mu] = PreFacP * ( XP[Mu] - QP[Mu] ); 
   };

  TDIData->UserFunction(X,  b,  3.0*PreFac,
                        XP, bP, 3.0*PreFacP,
                        TDIData->UserData, f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int yxIntegrand(unsigned ndim, const double *yx, void *params,
                unsigned fdim, double *fval)
{
  (void)ndim;
  TDIntegrandData *TDIData = (TDIntegrandData *)params;
  const double *y = yx + 0;
  const double *x = yx + 3;

  double u[NUMREGIONS][3];
  double ABCDEF[NUMREGIONS][6];
  double yJacobian[NUMREGIONS];
  GetSubregionData(y[0], y[1], y[2], u, ABCDEF, yJacobian);

  double *fThisRegion = TDIData->fBuffer;
  memset(fval, 0, fdim*sizeof(double));
  for(int d=0; d<NUMREGIONS; d++)
   { 
     if (fabs(yJacobian[d]) < 1.0e-8) 
      continue; 

     if ( RegionOnly!=-1 && RegionOnly!=d )
      continue;

     double Jacobian;
     GetxIntegrand(x, u[d], ABCDEF[d], TDIData, &Jacobian, fThisRegion);
     Jacobian *= yJacobian[d];

     if (fabs(Jacobian) < 1.0e-8) 
      continue;
     
     for(unsigned int nf=0; nf<fdim; nf++)
      fval[nf] += Jacobian*fThisRegion[nf];
   };

  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TetTetInt_TD(SWGVolume *OA, int ntA, int iQA,
                  SWGVolume *OB, int ntB, int iQB,
                  UserTTIntegrand UserIntegrand,
                  void *UserData, int fdim,
                  double *Result, double *Error,
                  double *fBuffer,
                  int MaxEvals, double RelTol)
{
  SWGTet *TA     = OA->Tets[ntA];
  double *QA     = OA->Vertices + 3*(TA->VI[ iQA ]);
  double PreFacA = (OA->Faces[TA->FI[iQA]]->Area) / (3.0 * TA->Volume);

  SWGTet *TB     = OB->Tets[ntB];
  double *QB     = OB->Vertices + 3*(TB->VI[ iQB ]);
  double PreFacB = (OB->Faces[TB->FI[iQB]]->Area) / (3.0 * TB->Volume);

  int OVIA[4], OVIB[4]; 
  double *V0A, *V1A, *V2A, *V3A, *V0B, *V1B, *V2B, *V3B;
  int ncv=CompareTets(OA, ntA, OB, ntB, OVIA, OVIB);
  V0A = V0B = OA->Vertices + 3*OVIA[0];
  V1A = V1B = OA->Vertices + 3*OVIA[1];
  V2A = V2B = OA->Vertices + 3*OVIA[2];
  V3A = V3B = OA->Vertices + 3*OVIA[3];
  if (ncv<4)
   V3B = OA->Vertices + 3*OVIB[3];
  if (ncv<3)
   V2B = OA->Vertices + 3*OVIB[2];
  if (ncv<2)
   V1B = OA->Vertices + 3*OVIB[1];
  if (ncv<1)
   V0B = OA->Vertices + 3*OVIB[0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TDIntegrandData MyData, *Data=&MyData;
  Data->UserFunction = UserIntegrand;
  Data->UserData     = UserData;
  Data->fBuffer      = fBuffer;
  
  Data->V0 = V0A;
  VecSub(V1A, V0A, Data->L1);
  VecSub(V2A, V1A, Data->L2);
  VecSub(V3A, V2A, Data->L3);
  Data->Q  = QA;
  Data->PreFac = PreFacA;

  Data->V0P = V0B;
  VecSub(V1B, V0B, Data->L1P);
  VecSub(V2B, V1B, Data->L2P);
  VecSub(V3B, V2B, Data->L3P);
  Data->QP  = QB;
  Data->PreFacP = PreFacB;

  double Lower[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Upper[6]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double Jacobian = 36.0 * TA->Volume * TB->Volume;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  hcubature(fdim, yxIntegrand, (void *)Data, 6, Lower, Upper,
	    MaxEvals, 0.0, RelTol, ERROR_INDIVIDUAL, Result, Error);

  for(int nf=0; nf<fdim; nf++)
   { Result[nf]*=Jacobian;
     Error[nf]*=Jacobian;
   };

}

} // namespace buff{
