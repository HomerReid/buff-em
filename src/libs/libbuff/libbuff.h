/* Copyright (C) 2005-2011 M. T. Homer Reid
 
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
 * libbuff.h -- header file for libbuff
 *
 * homer reid  -- 4/2014
 */

#ifndef LIBBUFF_H
#define LIBBUFF_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <complex>
#include <cmath>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libIncField.h>

#include "GTransformation.h"
#include "IHAIMatProp.h"

using namespace scuff;

namespace buff {

#ifndef ZVAC
#define ZVAC 376.73031346177
#endif

/***************************************************************/
/* SWGTet is a structure containing data on a single           */
/* tetrahedron in the meshed geometry.                         */
/***************************************************************/
typedef struct SWGTet
 { 
   int VI[4];                /* indices of vertices in Vertices array */
   int FI[4];                /* indices of faces in Faces array */
   double Centroid[3];
   double Volume;
   int Index;                /* index of this tet within SWGSurface (0..NumTets-1)*/

} SWGTet;

/***************************************************************/
/* SWGFace is a structure containing data on a single          */
/* triangular face in a meshed surface.                        */
/*                                                             */
/* This may be an *interior* face, in which case all fields in */
/* the structure are valid; or it may be an *exterior* face,   */
/* in which case the iQM, iMTet, and MIndex fields all have    */
/* the value -1.                                               */
/*                                                             */
/* note: after the following code snippet                      */
/*  double *R= O->Vertices + 3*F->IQP                          */
/* we have that R[0..2] are the cartesian coordinates of the   */
/* QP vertex in the basis function corresponding to face F.    */
/***************************************************************/
typedef struct SWGFace 
 { 
   int iV1, iV2, iV3, iQP, iQM;	/* indices of panel vertices (iV1<iV2<iV3) */
   double Centroid[3];          /* face centroid */
   double Area;                 /* area of face */
   double Radius;               /* radius of enclosing sphere */

   int iPTet;                   /* index of PTet within SWGVolume (0..NumTets-1)*/
   int iMTet;                   /* index of MTet within SWGVolume (0..NumTets-1)*/
   int PIndex;                  /* index of this face within PTet (0..3)*/
   int MIndex;                  /* index of this face within MTet (0..3)*/
   int Index;                   /* index of this face within SWGVolume (0..NumInteriorFaces-1)*/

   SWGFace *Next;               /* pointer to next face in linked list */

} SWGFace;

/***************************************************************/
/***************************************************************/
/***************************************************************/
class SWGVolume
 { 
   /*--------------------------------------------------------------*/ 
   /*- class methods ----------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
  public:  

   /*-------------------------------------------------------------------*/
   /* main constructor entry point: construct from 'OBJECT...ENDOBJECT' */
   /* or SURFACE...ENDSURFACE section in a .scuffgeo file               */
   /*-------------------------------------------------------------------*/
   SWGVolume(char *MeshFileName, char *Label=0, char *MatFileName=0,
             GTransformation *OTGT=0);

   /* destructor */
   ~SWGVolume();

   /*-------------------------------------------------------------------*/
   /*- visualization ---------------------------------------------------*/
   /*-------------------------------------------------------------------*/
   void DrawBF(int nf, FILE *f);
   void DrawFace(int nf, FILE *f);
   void DrawTet(int nt, FILE *f);

   /*-------------------------------------------------------------------*/
   /*- apply geometrical transformations -------------------------------*/
   /*-------------------------------------------------------------------*/
   void Transform(const GTransformation *DeltaGT);
   void Transform(const char *format,...);
   void UnTransform();

//  private:

   /*--------------------------------------------------------------*/
   /*- private data fields  ---------------------------------------*/
   /*--------------------------------------------------------------*/
   int NumVertices;                /* number of vertices in mesh  */
   double *Vertices;               /* Vertices[3*n,3*n+1,3*n+2]=nth vertex coords */

   int NumTets;                    /* number of tetrahedra */
   SWGTet **Tets;                  /* array of pointers to tetrahedra */

   int NumInteriorFaces;           /* number of interior triangular faces */
   int NumTotalFaces;              /* number of total triangular faces */
   SWGFace **Faces;                /* array of pointers to interior faces */
   SWGFace **ExteriorFaces;        /* array of pointers to exterior faces */

   char *MeshFileName;             /* saved name of mesh file */
   char *Label;                    /* unique label identifying the volume*/

   char *MatFileName;              /* saved name of material property file */
   IHAIMatProp *MP;                /* material property */

   /* GT encodes any transformation that has been carried out since */
   /* the object was read from its mesh file (not including a       */
   /* possible one-time GTransformation that may have been specified*/
   /* in the .scuffgeo file when the surface was first created.)    */
   GTransformation *GT;
   GTransformation *OTGT;

   // the following fields are used to pass some data items up to the
   // higher-level routine that calls the SWGVolume constructor
   char *ErrMsg;  /* used to indicate error to calling routine */

   /*--------------------------------------------------------------*/ 
   /*- private class methods --------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /* constructor subroutines */
   void InitFaceList();
   void ReadGMSHFile(FILE *MeshFile);
   void InitSWGFace(SWGFace *F);

 }; // class SWGVolume

/***************************************************************/
/***************************************************************/
/***************************************************************/
class SWGGeometry
 { 
   /*--------------------------------------------------------------*/ 
   /*- class methods ----------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
  public:  

   // constructor 
   SWGGeometry(const char *GeoFileName);

   // destructor
   ~SWGGeometry();

   // geometric transforms 
   void Transform(GTComplex *GTC);
   void UnTransform();
   char *CheckGTCList(GTComplex **GTCList, int NumGTCs);

   // visualization
   void WritePPMesh(const char *FileName, const char *Tag);
   void PlotPermittivity(const char *FileName, const char *Tag);
   void PlotCurrentDistribution(const char *FileName, HVector *J,
                                const char *Tag, ...);

   // scattering API
   HMatrix *AllocateVIEMatrix(bool PureImagFreq=false);
   HMatrix *AssembleVIEMatrix(cdouble Omega, HMatrix *M);
   HVector *AllocateRHSVector();
   HVector *AssembleRHSVector(cdouble Omega, IncField *IF, HVector *V);
   void GetFields(IncField *IF, HVector *J, cdouble Omega, double *X, cdouble *EH);
   HMatrix *GetFields(IncField *IF, HVector *J, cdouble Omega,
                      HMatrix *XMatrix, HMatrix *FMatrix=NULL);
   HMatrix *GetPFT(IncField *IF, HVector *JVector,
                   cdouble Omega, HMatrix *PFTMatrix=0, bool *NeedQuantity=0);

   // compute individual matrix blocks
   void AssembleGBlock(int noa, int nob, cdouble Omega,
                       HMatrix *G,
                       int RowOffset=0, int ColOffset=0);
   void AssembleVInvBlock(int no, cdouble Omega,
                          SMatrix *VInv, SMatrix *ImEps,
                          HMatrix *TInv=0, int Offset=0);

   // miscellaneous routines
   SWGVolume *GetObjectByLabel(const char *Label, int *pno=0);

   // directories within which to search for mesh files
   static int NumMeshDirs;
   static char **MeshDirs;
   static int MaxTDEvals;
   static int TDNCVThreshold;

//  private:

   /*--------------------------------------------------------------*/
   /*- private data fields  ---------------------------------------*/
   /*--------------------------------------------------------------*/
   int NumObjects;
   SWGVolume **Objects;

   int TotalBFs, *BFIndexOffset;

   int *Mate;

   char *GeoFileName;

 }; // class SWGGeometry

/***************************************************************/
/* non-class utility methods ***********************************/
/***************************************************************/
SWGTet *NewSWGTet(double *Vertices, int iV1, int iV2, int iV3, int iV4);

int CompareBFs(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB,
               double *rRel=0);

int CompareTets(SWGVolume *OA, int ntA, SWGVolume *OB, int ntB,
                int *OVIA=0, int *OVIB=0);

int GetOverlapElements(SWGVolume *O, int nfA,
                       int Indices[7], double Entries[7]);

void GetDQMoments(SWGVolume *O, int nf, double J[3], double Q[3][3],
                  bool NeedQ=true);

/***************************************************************/
/* routine to compute matrix elements of the dyadic GF and its */
/* derivatives                                                 */
/***************************************************************/
cdouble GetGMatrixElement(SWGVolume *VA, int nfA,
                          SWGVolume *VB, int nfB,
                          cdouble Omega,
                          bool *NeedDerivatives=0,
                          cdouble *dG=0,
                          int rPower=-10, bool ForceBF=false);

/***************************************************************/
/* routines for integrating over tetrahedra and faces          */
/***************************************************************/
typedef void (*UserTIntegrand)(double *x, double *b, double Divb,
                               void *UserData, double *I);

typedef void (*UserTTIntegrand)(double *xA, double *bA, double DivbA,
                                double *xB, double *bB, double DivbB,
                                void *UserData, double *I);

typedef void (*UserFIntegrand)(double *x, double *b, double Divb, double *nHat,
                               void *UserData, double *I);

typedef void (*UserFFIntegrand)(double *xA, double *bA, double DivbA, double *nHatA,
                                double *xB, double *bB, double DivbB, double *nHatB,
                                void *UserData, double *I);

void TetInt(SWGVolume *V, int nt, int iQ, double Sign,
            UserTIntegrand Integrand, void *UserData,
            int fdim, double *Result, double *Error,
            int NumPts, int MaxEvals, double RelTol);

void BFInt(SWGVolume *V, int nf,
           UserTIntegrand Integrand, void *UserData,
           int fdim, double *Result, double *Error,
           int NumPts, int MaxEvals, double RelTol);

void TetTetInt(SWGVolume *VA, int ntA, int iQA, double SignA,
               SWGVolume *VB, int ntB, int iQB, double SignB,
               UserTTIntegrand Integrand, void *UserData,
               int fdim, double *Result, double *Error,
               int NumPts, int MaxEvals, double RelTol);

void TetTetInt_TD(SWGVolume *VA, int ntA, int iQA,
                  SWGVolume *VB, int ntB, int iQB,
                  UserTTIntegrand UserIntegrand,
                  void *UserData, int fdim,
                  double *Result, double *Error, double *Buffer,
                  int MaxEvals, double RelTol);

void TetTetInt_SI(SWGVolume *VA, int ntA, int iQa,
                  SWGVolume *VB, int ntB, int iQb,
                  UserFFIntegrand Integrand, void *UserData,
                  int fdim, double *Result, double *Error,
                  int Order, int MaxEvals, double RelTol);

void BFBFInt(SWGVolume *VA, int nfA,
             SWGVolume *VB, int nfB,
             UserTTIntegrand Integrand, void *UserData,
             int fdim, double *Result, double *Error,
             int NumPts, int MaxEvals, double RelTol);

void FaceInt(SWGVolume *V, int nt, int nf, int nfBF, double Sign,
             UserFIntegrand Integrand, void *UserData,
             int fdim, double *Result, double *Error,
             int Order, int MaxEvals, double RelTol);

void FaceFaceInt(SWGVolume *VA, int ntA, int nfA, int nfBFA, double SignA,
                 SWGVolume *VB, int ntB, int nfB, int nfBFB, double SignB,
                 UserFFIntegrand Integrand, void *UserData,
                 int fdim, double *Result, double *Error,
                 int Order, int MaxEvals, double RelTol);

double *GetTetCR(int NumPts);

} // namespace buff 

#endif // #ifndef LIBBUFF_H
