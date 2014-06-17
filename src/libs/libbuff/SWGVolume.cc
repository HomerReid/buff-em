/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * SWGVolume.cc -- implementation of some methods in the SWGVolume
 *               -- class (this is the class formerly known as RWGObject)
 *
 * homer reid    -- 3/2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "cmatheval.h"
#include "libbuff.h"
#include "GTransformation.h"

using namespace scuff;

namespace buff {

#define MAXSTR 1000
#define MAXTOK 50  

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
SWGTet *NewSWGTet(double *Vertices, int iV1, int iV2, int iV3, int iV4)
{ 
  SWGTet *T = (SWGTet *)mallocEC(sizeof *T);
  T->VI[0] = iV1;
  T->VI[1] = iV2;
  T->VI[2] = iV3;
  T->VI[3] = iV4;

  // T->FI[n] is the index within the Faces array of the
  // face opposite vertex #n; these are not known until
  // after InitFaceList().
  T->FI[0] = T->FI[1] = T->FI[2] = T->FI[3] = -1;

  // compute centroid 
  double *V1 = Vertices + 3*iV1; 
  double *V2 = Vertices + 3*iV2; 
  double *V3 = Vertices + 3*iV3; 
  double *V4 = Vertices + 3*iV4; 
  for(int i=0; i<3; i++)
   T->Centroid[i] = 0.25*( V1[i] + V2[i] + V3[i] + V4[i] );

  // volume = (1/6) A \cdot (B x C)
  double A[3], B[3], C[3];
  for(int i=0; i<3; i++)
   { A[i] = V2[i] - V1[i];
     B[i] = V3[i] - V1[i];
     C[i] = V4[i] - V1[i];
   };
  T->Volume = (  A[0] * (B[1]*C[2] - B[2]*C[1])
                +A[1] * (B[2]*C[0] - B[0]*C[2])
                +A[2] * (B[0]*C[1] - B[1]*C[0])
              ) / 6.0;

  return T;

}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
SWGVolume::SWGVolume(char *pMeshFileName,
                     char *pLabel,
                     char *pMatFileName,
                     GTransformation *pOTGT)
{
  ErrMsg=0;

  /*------------------------------------------------------------*/
  /*- try to open the mesh file. we look in several places:     */
  /*- (a) the current working directory                         */
  /*- (b) any directories that may have been specified by       */
  /*-     MESHPATH statements in .scuffgeo files or via the     */
  /*-     SCUFF_MESH_PATH environment variable                  */
  /*------------------------------------------------------------*/
  MeshFileName=strdup(pMeshFileName);
  FILE *MeshFile=fopen(MeshFileName,"r");
  if (!MeshFile)
   { for(int nmd=0; MeshFile==0 && nmd<RWGGeometry::NumMeshDirs; nmd++)
      { MeshFile=vfopen("%s/%s","r",RWGGeometry::MeshDirs[nmd],MeshFileName);
        if (MeshFile) 
         Log("Found mesh file %s/%s",RWGGeometry::MeshDirs[nmd],MeshFileName);
      };
   };
  if (!MeshFile)
   ErrExit("could not open file %s",MeshFileName);
   
  /*------------------------------------------------------------*/
  /*- initialize simple fields ---------------------------------*/
  /*------------------------------------------------------------*/
  NumVertices=NumTets=NumInteriorFaces=NumTotalFaces=0;
  Vertices=0;
  Tets=0;
  Faces=0;
  if (pLabel==0)
   Label=strdup(MeshFileName);
  else
   Label=strdup(pLabel);

  if (pMatFileName==0)
   { MatFileName=0;
     MP=0;
   }
  else
   { MatFileName=strdup(pMatFileName);
     MP=new IHAIMatProp(MatFileName);
     if (MP->ErrMsg)
      ErrExit(MP->ErrMsg);
   };

  /*------------------------------------------------------------*/
  /*- note: the 'OTGT' parameter to this function is distinct   */
  /*- from the 'GT' field inside the class body. the former is  */
  /*- an optional 'One-Time Geometrical Transformation' to be   */
  /*- applied to the object once at its creation. the latter    */
  /*- is designed to store a subsequent transformation that may */
  /*- be applied to the surface, and is initialized to zero.    */
  /*-                                                           */
  /*- Note: in contrast to the implementation of RWGSurface in  */
  /*- libscuff, here we need to hold on to the initial OTGT *as */
  /*- well as* the running GT. The reason is that, when we      */
  /*- untransform the coordinates of points inside the volume   */
  /*- for the purposes of evaluating the position-dependent     */
  /*- permittivity tensor, we need to transform all the way     */
  /*- back to the untransformed object as described in the mesh */
  /*- file, not just to the partially transformed object as     */
  /*- described in the .buffgeo file.                           */
  /*------------------------------------------------------------*/
  GT=0;
  OTGT=pOTGT;

  /*------------------------------------------------------------*/
  /*- Switch off based on the file type to read the mesh file:  */
  /*-  1. file extension=.msh    --> ReadGMSHFile              -*/
  /*-  2. otherwise error (for now)                            -*/
  /*------------------------------------------------------------*/
  char *p=GetFileExtension(MeshFileName);
  if (!p)
   ErrExit("file %s: invalid extension",MeshFileName);
  else if (!StrCaseCmp(p,"msh"))
   ReadGMSHFile(MeshFile,0);
  else
   ErrExit("file %s: unknown extension %s",MeshFileName,p);

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (NumTets==0)
   ErrExit("file %s: no tetrahedra found",MeshFileName);

  /*------------------------------------------------------------*/
  /* gather necessary edge connectivity info. this is           */
  /* complicated enough to warrant its own separate routine.    */
  /*------------------------------------------------------------*/
  InitFaceList();

  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  /*------------------------------------------------------------*/
  if (OTGT)
   { OTGT->Apply(Vertices, NumVertices);
     for(int nf=0; nf<NumTotalFaces; nf++)
      OTGT->Apply(Faces[nf]->Centroid, 1);
     for(int nt=0; nt<NumTets; nt++)
      OTGT->Apply(Tets[nt]->Centroid, 1);
   };

} 

/***************************************************************/
/* SWGVolume destructor.                                       */
/***************************************************************/
SWGVolume::~SWGVolume()
{ 
  free(Vertices);

  for(int nf=0; nf<NumTotalFaces; nf++)
   free(Faces[nf]);
  free(Faces);

  for(int nt=0; nt<NumTets; nt++)
   free(Tets[nt]);
  free(Tets);

  if (MeshFileName) free(MeshFileName);
  if (Label) free(Label);
  if (GT) delete GT;
  if (ErrMsg) free(ErrMsg);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGVolume::Transform(const GTransformation *DeltaGT)
{ 
  /***************************************************************/
  /*- first apply the transformation to all points whose         */
  /*- coordinates we store inside the SWGVolume structure:       */
  /*- vertices, face centroids, and tet centroids.               */
  /***************************************************************/
  /* vertices */
  DeltaGT->Apply(Vertices, NumVertices);

  /* face centroids */
  for(int nf=0; nf<NumTotalFaces; nf++)
   DeltaGT->Apply(Faces[nf]->Centroid, 1);

  /* tet centroids */
  for(int nt=0; nt<NumTets; nt++)
   DeltaGT->Apply(Tets[nt]->Centroid, 1);

  /***************************************************************/
  /* update the internally stored GTransformation ****************/
  /***************************************************************/
  if (!GT)
    GT = new GTransformation(DeltaGT);
  else
    GT->Transform(DeltaGT);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SWGVolume::Transform(const char *format,...)
{
  va_list ap;
  char buffer[MAXSTR];
  va_start(ap,format);
  vsnprintfEC(buffer,MAXSTR,format,ap);

  GTransformation MyGT(buffer, &ErrMsg);
  if (ErrMsg)
   ErrExit(ErrMsg);
  Transform(&MyGT);
}

/***************************************************************/
/* undo the internally-stored GTransformation                  */
/***************************************************************/
void SWGVolume::UnTransform()
{
  if (!GT)
   return;
 
  /***************************************************************/
  /* untransform vertices                                        */
  /***************************************************************/
  GT->UnApply(Vertices, NumVertices);

  /***************************************************************/
  /* untransform face / tet centroids                            */
  /***************************************************************/
  for(int nf=0; nf<NumTotalFaces; nf++)
   GT->UnApply(Faces[nf]->Centroid, 1);
  for(int nt=0; nt<NumTets; nt++)
   GT->UnApply(Tets[nt]->Centroid, 1);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GT->Reset();

}

/***************************************************************/
/* compare two SWG basis functions by computing the following  */
/* characteristics:                                            */
/*  (a) number of common vertices                              */
/*  (b) relative distance                                      */
/***************************************************************/
int CompareBFs(SWGVolume *OA, int nfA, SWGVolume *OB, int nfB,
               double *rRel)
{ 

  SWGFace *FA = OA->Faces[nfA];
  SWGFace *FB = OB->Faces[nfB];

  if (rRel)
   { *rRel = VecDistance(FA->Centroid, FB->Centroid) 
              / fmax(FA->Radius, FB->Radius); 
   };

  if (OA!=OB)
   return 0;
 
  int ncv=0;

  if (FA->iQP == FB->iQP ) ncv++;
  if (FA->iQP == FB->iV1 ) ncv++;
  if (FA->iQP == FB->iV2 ) ncv++;
  if (FA->iQP == FB->iV3 ) ncv++;
  if (FA->iQP == FB->iQM ) ncv++;

  if (FA->iV1 == FB->iQP ) ncv++;
  if (FA->iV1 == FB->iV1 ) ncv++;
  if (FA->iV1 == FB->iV2 ) ncv++;
  if (FA->iV1 == FB->iV3 ) ncv++;
  if (FA->iV1 == FB->iQM ) ncv++;

  if (FA->iV2 == FB->iQP ) ncv++;
  if (FA->iV2 == FB->iV1 ) ncv++;
  if (FA->iV2 == FB->iV2 ) ncv++;
  if (FA->iV2 == FB->iV3 ) ncv++;
  if (FA->iV2 == FB->iQM ) ncv++;

  if (FA->iV3 == FB->iQP ) ncv++;
  if (FA->iV3 == FB->iV1 ) ncv++;
  if (FA->iV3 == FB->iV2 ) ncv++;
  if (FA->iV3 == FB->iV3 ) ncv++;
  if (FA->iV3 == FB->iQM ) ncv++;

  if (FA->iQM == FB->iQP ) ncv++;
  if (FA->iQM == FB->iV1 ) ncv++;
  if (FA->iQM == FB->iV2 ) ncv++;
  if (FA->iQM == FB->iV3 ) ncv++;
  if (FA->iQM == FB->iQM ) ncv++;

  return ncv;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int GetOverlapElements(SWGVolume *O, int nfA,
                       int Indices[7], double Entries[7])
{ 
  int NNZ=0;
  
  SWGFace *FA  = O->Faces[nfA];
  double *QP   = O->Vertices + 3*(FA->iQP);
  double *V1   = O->Vertices + 3*(FA->iV1);
  double *V2   = O->Vertices + 3*(FA->iV2);
  double *V3   = O->Vertices + 3*(FA->iV3);
  double *QM   = O->Vertices + 3*(FA->iQM);
  SWGTet *PTet = O->Tets[ FA->iPTet ];
  SWGTet *MTet = O->Tets[ FA->iMTet ];

  // self-overlap term is always the first entry in the output arrays
  Indices[0] = nfA;
  Entries[0] = 0.0;
  NNZ=1;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int SignA = 1; SignA >= -1; SignA-=2)
   { 
     SWGTet *T  = (SignA==1) ? PTet : MTet;
     double *QA = (SignA==1) ? QP : QM;

     double L1[3], L2[3], L3[3], XmQ[3];
     VecSub(V1, QA, L1);
     VecSub(V2, QA, L2);
     VecSub(V3, QA, L3);
     VecSub(FA->Centroid, QA, XmQ);
     double XmQ2=VecNorm2(XmQ);
     double L2Sum=VecNorm2(L1) + VecNorm2(L2) + VecNorm2(L3);

     // diagonal contribution
     Entries[0] += (FA->Area)*(FA->Area)*(XmQ2/20.0 + L2Sum/180.0) / T->Volume;

     // off-diagonal contributions
     for(int iF=0; iF<4; iF++)
      {
        int nfB = T->FI[iF];
        if ( nfB == nfA ) continue;
        if ( nfB >= O->NumInteriorFaces ) continue;

        SWGFace *FB = O->Faces[ nfB ];
        double SignB, *QB;
        if ( FB->iPTet == T->Index )
         { SignB = 1.0;
           QB = O->Vertices + 3*FB->iQP;
         }
        else
         { SignB=-1.0;
           QB = O->Vertices + 3*FB->iQM;
         };

        double PreFac = 2.0*SignA*SignB*(FA->Area)*(FB->Area)
                       / (3.0*T->Volume);
        double XmQdDQ=   XmQ[0]*(QA[0]-QB[0])
                       + XmQ[1]*(QA[1]-QB[1])
                       + XmQ[2]*(QA[2]-QB[2]);

        Indices[NNZ] = nfB;
        Entries[NNZ] = PreFac * (  3.0*XmQ2/40.0
                                  +L2Sum/120.0
                                  +XmQdDQ/8.0
                                );
        NNZ++;

      }; // for(int iF=0; iF<4; iF++)

   }; // for(int SignA = 1; SignA >=-1; SignA-=2)

  return NNZ;

}

} // namespace buff
