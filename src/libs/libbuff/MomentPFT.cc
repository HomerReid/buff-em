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
 * MomentPFT.cc -- BUFF-EM code module to compute lowest-order
 *                 multipole contributions to power, force, and torque
 *
 * homer reid   -- 8/2015
 *
 */

#include <libscuff.h>
#include <libbuff.h>

#define II cdouble(0.0,1.0)

namespace buff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Get1BFMoments(SWGVolume *O, int nf,
                   double JMu[3], double JMuNu[3][3])
{
  SWGFace *F  = O->Faces[nf];
  double  A   = F->Area;
  double *QP  = O->Vertices + 3*(F->iQP);
  double *QM  = O->Vertices + 3*(F->iQM);
  
  for(int Mu=0; Mu<3; Mu++)
   JMu[Mu] = A*(QM[Mu] - QP[Mu]) / 4.0;

// note: a more compact formula for the following would be nice
  double *V1  = O->Vertices + 3*(F->iV1);
  double *V2  = O->Vertices + 3*(F->iV2);
  double *V3  = O->Vertices + 3*(F->iV3);
  double L1P[3], L2P[3], L3P[3], L1M[3], L2M[3], L3M[3];
  VecSub(V1, QP, L1P);
  VecSub(V2, QP, L2P);
  VecSub(V3, QP, L3P);
  VecSub(V1, QM, L1M);
  VecSub(V2, QM, L2M);
  VecSub(V3, QM, L3M);

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    JMuNu[Mu][Nu] 
     = + A*(  L1P[Mu]*L1P[Nu] - L1M[Mu]*L1M[Nu]
            + L2P[Mu]*L2P[Nu] - L2M[Mu]*L2M[Nu]
            + L3P[Mu]*L3P[Nu] - L3M[Mu]*L3M[Nu]
           ) / 30.0
       + A*(  L1P[Mu]*L2P[Nu] - L1M[Mu]*L2M[Nu]
             +L1P[Mu]*L3P[Nu] - L1M[Mu]*L3M[Nu]
             +L2P[Mu]*L3P[Nu] - L2M[Mu]*L3M[Nu]
             +L1P[Nu]*L2P[Mu] - L1M[Nu]*L2M[Mu]
             +L1P[Nu]*L3P[Mu] - L1M[Nu]*L3M[Mu]
             +L2P[Nu]*L3P[Mu] - L2M[Nu]*L3M[Mu]
           ) /60.0
       + A*(  L1P[Mu]*QP[Nu]  - L1M[Mu]*QM[Nu]
            + L2P[Mu]*QP[Nu]  - L2M[Mu]*QM[Nu]
            + L3P[Mu]*QP[Nu]  - L3M[Mu]*QM[Nu]
           ) /12.0;
}

/***************************************************************/
/* if Workspace is non-null it must point to a buffer with     */
/* space for at least 12*N cdoubles, where N is the number of  */
/* basis functions on object #no.                              */
/***************************************************************/
void GetNEQMoments(SWGGeometry *G, int no, HMatrix *Rytov,
                   cdouble MMuNu[3][3], cdouble MMuNuRho[3][3][3],
                   cdouble *Workspace=0)
{
  SWGVolume *O      = G->Objects[no];
  int N             = O->NumInteriorFaces;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  static int SavedN=0;
  static cdouble *SavedWorkspace=0;
  if ( Workspace==0 )
   { if (SavedN==N)
      Workspace = SavedWorkspace;
     else
      { if (SavedWorkspace) 
         free(SavedWorkspace);
        SavedN=N;
        SavedWorkspace = (cdouble *)mallocEC(12*N*sizeof(cdouble));
        Workspace = SavedWorkspace;
      };
   };
    
  int nVec=0;
  HVector *vMu[3], *vMuNu[3][3];
  for(int Mu=0; Mu<3; Mu++, nVec++)
   vMu[Mu] = new HVector(N, LHM_COMPLEX, Workspace + nVec*N);
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++, nVec++)
    vMuNu[Mu][Nu] = new HVector(N, LHM_COMPLEX, Workspace + nVec*N);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("GNM Getting 1BF moments...");
  for(int n=0; n<N; n++)
   { 
     double JMu[3], JMuNu[3][3];
     Get1BFMoments(O, n, JMu, JMuNu);

     for(int Mu=0; Mu<3; Mu++)
      vMu[Mu]->SetEntry(n, JMu[Mu]);

     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       vMuNu[Mu][Nu]->SetEntry(n, JMuNu[Mu][Nu]);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("GNM doing linear algebra...");

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    MMuNu[Mu][Nu]=Rytov->BilinearProduct(vMu[Mu], vMu[Nu]);

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    for(int Rho=0; Rho<3; Rho++)
     MMuNuRho[Mu][Nu][Rho]=Rytov->BilinearProduct(vMuNu[Mu][Rho], vMu[Nu]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int Mu=0; Mu<3; Mu++)
   delete vMu[Mu];
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    delete vMuNu[Mu][Nu];
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetMoments(SWGGeometry *G, int no, cdouble Omega, 
                HVector *JVector,
                cdouble p[3], cdouble m[3], cdouble Qp[3])
{
  p[0]=p[1]=p[2]=m[0]=m[1]=m[2]=0.0;
  cdouble Q[3][3]={{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  cdouble iw=II*Omega;
  SWGVolume *O=G->Objects[no];
  int Offset = G->BFIndexOffset[no];
  for(int n=0; n<O->NumInteriorFaces; n++)
   { 
     cdouble JAlpha = JVector->GetEntry(Offset + n);

     double JMu[3], JMuNu[3][3];
     Get1BFMoments(O, n, JMu, JMuNu);

     cdouble JMMTrace=JAlpha*(JMuNu[0][0]+JMuNu[1][1]+JMuNu[2][2]);
     for(int Mu=0; Mu<3; Mu++)
      { 
        p[Mu] -= JAlpha*JMu[Mu] / iw;

        int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
        m[Mu] -= 0.5*JAlpha*(JMuNu[MP1][MP2]-JMuNu[MP2][MP1]);

        Q[Mu][Mu] -= JAlpha*(6.0*JMuNu[Mu][Mu] - 2.0*JMMTrace)/iw;
        for(int Nu=Mu+1; Nu<3; Nu++)
         Q[Mu][Nu] -= 3.0*JAlpha*(JMuNu[Mu][Nu] + JMuNu[Nu][Mu])/iw;
      };
     Q[1][0]=Q[0][1];
     Q[2][0]=Q[0][2];
     Q[2][1]=Q[1][2];
   };

  Qp[0]=Qp[1]=Qp[2]=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    Qp[Mu] += Q[Mu][Nu] * p[Nu];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DoPrincipalAxisDecomposition(cdouble MMuNu[3][3], 
                                  cdouble MMuNuRho[3][3][3], 
                                  cdouble Omega,
                                  cdouble p[3][3], cdouble m[3][3])
{
  static HMatrix *ppMatrix = new HMatrix(3,3,LHM_COMPLEX);
  static HMatrix *mpMatrix = new HMatrix(3,3,LHM_COMPLEX);
  static HMatrix *U        = new HMatrix(3,3,LHM_COMPLEX);
  static HMatrix *VT       = new HMatrix(3,3,LHM_COMPLEX);
  static HVector *Lambda   = new HVector(3,LHM_REAL);
  static HVector *Sigma    = new HVector(3,LHM_REAL);

  double w2=real(Omega)*real(Omega);
  for(int Mu=0; Mu<3; Mu++)
   { ppMatrix->SetEntry(Mu, Mu, real(MMuNu[Mu][Mu])/w2);
     for(int Nu=Mu+1; Nu<3; Nu++)
      { ppMatrix->SetEntry(Mu, Nu, MMuNu[Mu][Nu]       / w2 );
        ppMatrix->SetEntry(Nu, Mu, conj(MMuNu[Mu][Nu]) / w2 );
      };
   };

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { int NP1=(Nu+1)%3, NP2=(Nu+2)%3;
      mpMatrix->SetEntry(Mu,Nu, II*Omega*(  MMuNuRho[Mu][NP1][NP2]
                                          - MMuNuRho[Mu][NP2][NP1]
                                         )
                        );
    };

  ppMatrix->Eig(Lambda, U);
  double MaxEig=abs(Lambda->GetEntry(0));
  MaxEig=fmax(MaxEig, abs(Lambda->GetEntry(1)));
  MaxEig=fmax(MaxEig, abs(Lambda->GetEntry(2)));

  for(int Mu=0.0; Mu<3.0; Mu++)
   if( abs(Lambda->GetEntry(Mu)) < 1.0e-8*MaxEig )
    Lambda->SetEntry(Mu, 0.0);

  for(int a=0; a<3; a++)
   for(int Mu=0; Mu<3; Mu++)
    p[a][Mu] = sqrt(Lambda->GetEntry(a))*U->GetEntry(Mu,a);

  mpMatrix->SVD(Sigma, U, VT);
  for(int a=0; a<3; a++)
   { 
     cdouble DotProd = conj(VT->GetEntry(a,0)) * p[a][0]
                      +conj(VT->GetEntry(a,1)) * p[a][1]
                      +conj(VT->GetEntry(a,2)) * p[a][2];
     cdouble LambdaA=Lambda->GetEntry(a);
     cdouble ScaleFactor
      = (LambdaA==0.0) ? 0.0 : Sigma->GetEntry(a) * DotProd/ LambdaA;
     m[a][0] = ScaleFactor*U->GetEntry(0,a);
     m[a][1] = ScaleFactor*U->GetEntry(1,a);
     m[a][2] = ScaleFactor*U->GetEntry(2,a);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetMomentPFT(SWGGeometry *G, int no, cdouble Omega,
                  HVector *JVector, HMatrix *Rytov,
                  HMatrix *PFTMatrix, bool WritePPFile, double QPF[3])
{
  // p[a][Mu] = Muth component of ath principal p-vector
  // m[a][Mu] = Muth component of ath principal m-vector
  cdouble p[3][3];
  cdouble m[3][3];
  cdouble Qp[3];  // Q*p, only used in the deterministic case
  int NumMoments; // = 1 for deterministic, 3 for neq

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (JVector)
   { NumMoments=1;
     GetMoments(G, no, Omega, JVector, p[0], m[0], Qp);
   }
  else
   { 
     NumMoments=3;
     Log("GMP getting NEQ moments...");
     cdouble MMuNu[3][3], MMuNuRho[3][3][3];
     GetNEQMoments(G, no, Rytov, MMuNu, MMuNuRho);
     DoPrincipalAxisDecomposition(MMuNu, MMuNuRho, Omega, p, m);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double k3=real(Omega)*real(Omega)*real(Omega);
  double k4=real(Omega)*k3;
  double k5=real(Omega)*k4;
  double PPF  = k4*ZVAC/(12.0*M_PI);
  double FPF1 = -TENTHIRDS*k4*ZVAC/(12.0*M_PI);
  double FPF2 = -TENTHIRDS*k5*ZVAC/(120.0*M_PI);
  double TPF  = -TENTHIRDS*k3*ZVAC/(6.0*M_PI);

  PFTMatrix->Zero();
  for(int nm=0; nm<NumMoments; nm++)
   for(int Mu=0; Mu<3; Mu++)
    { 
      PFTMatrix->AddEntry(no, PFT_PABS, 
                              PPF*real( conj(p[nm][Mu])*p[nm][Mu] )
                         );

      int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
      PFTMatrix->AddEntry(no, PFT_XFORCE + Mu,
                              FPF1*real( conj(m[nm][MP1])*p[nm][MP2]
                                        -conj(m[nm][MP2])*p[nm][MP1] ) 
                         );

      PFTMatrix->AddEntry(no, PFT_XTORQUE + Mu,
                              TPF*( real(p[nm][MP1]) * imag(p[nm][MP2])
                                   -real(p[nm][MP2]) * imag(p[nm][MP1])
                                  )
                         );
    };

  if (JVector && QPF)
   { QPF[0] = FPF2*imag(Qp[0]);
     QPF[1] = FPF2*imag(Qp[1]);
     QPF[2] = FPF2*imag(Qp[2]);
   };

  if (WritePPFile)
   {
     double Radius=0.0;
     SWGVolume *O=G->Objects[no];
     for(int nv=0; nv<O->NumVertices; nv++)
      Radius=fmax(Radius, VecNorm(O->Vertices + 3*nv));
  
     FILE *f=vfopen("%s.Moments.pp","a",O->Label);

     for(int nm=0; nm<NumMoments; nm++)
      { 
        cdouble *pn=p[nm], *mn=m[nm];

        double pNorm = sqrt( norm(pn[0]) + norm(pn[1]) + norm(pn[2]) );
        double mNorm = sqrt( norm(mn[0]) + norm(mn[1]) + norm(mn[2]) );
        double pScaleFac = 1.5*Radius/pNorm;
        double mScaleFac = 1.5*Radius/mNorm;

        fprintf(f,"View \"Real p%i_%g\" {\n",nm+1,real(Omega));
        fprintf(f,"VP(0,0,0) {%e,%e,%e};\n",
                   pScaleFac*real(pn[0]),
                   pScaleFac*real(pn[1]),
                   pScaleFac*real(pn[2]));
        fprintf(f,"};\n");

        fprintf(f,"View \"Imag p%i_%g\" {\n",nm+1,real(Omega));
        fprintf(f,"VP(0,0,0) {%e,%e,%e};\n",
                   pScaleFac*imag(pn[0]),
                   pScaleFac*imag(pn[1]),
                   pScaleFac*imag(pn[2]));
        fprintf(f,"};\n");

        fprintf(f,"View \"Real m%i_%g\" {\n",nm+1,real(Omega));
        fprintf(f,"VP(0,0,0) {%e,%e,%e};\n",
                   mScaleFac*real(mn[0]),
                   mScaleFac*real(mn[1]),
                   mScaleFac*real(mn[2]));
        fprintf(f,"};\n");

        fprintf(f,"View \"Imag m%i_%g\" {\n",nm+1,real(Omega));
        fprintf(f,"VP(0,0,0) {%e,%e,%e};\n",
                   mScaleFac*imag(mn[0]),
                   mScaleFac*imag(mn[1]),
                   mScaleFac*imag(mn[2]));
        fprintf(f,"};\n");

      }; // for(int nm=0; nm<NumMoments; nm++)
    fclose(f);

   };

}

} // namespace buff
