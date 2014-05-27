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
 * InitFaceList.cc
 *
 * homer reid   -- 4/2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;
namespace buff{

/***************************************************************/
/* this routine just computes some geometric data for a given  */
/* triangular face.                                            */
/***************************************************************/
void SWGVolume::InitSWGFace(SWGFace *F)
{
  double *V1 = Vertices + 3*(F->iV1);
  double *V2 = Vertices + 3*(F->iV2);
  double *V3 = Vertices + 3*(F->iV3);

  for(int i=0; i<3; i++)
   F->Centroid[i]=(V1[i] + V2[i] + V3[i]) / 3.0;

  double A[3], B[3], AxB[3];
  A[0] = V2[0]-V1[0];
  A[1] = V2[1]-V1[1];
  A[2] = V2[2]-V1[2];
  B[0] = V3[0]-V1[0];
  B[1] = V3[1]-V1[1];
  B[2] = V3[2]-V1[2];
  AxB[0] = A[1]*B[2] - A[2]*B[1];
  AxB[1] = A[2]*B[0] - A[0]*B[2];
  AxB[2] = A[0]*B[1] - A[1]*B[0];
  F->Area = 0.5*sqrt( AxB[0]*AxB[0] + AxB[1]*AxB[1] + AxB[2]*AxB[2] );

}

/***************************************************************/
/* InitFaceList: After reading in a set of tetrahedra, we call */
/* this function to extract all necessary information          */
/* regarding interior faces, exterior faces, etc.              */
/***************************************************************/
void SWGVolume::InitFaceList()
{ 
  /***************************************************************/
  /***************************************************************/
  /*- allocate temporary storage arrays  *************************/
  /***************************************************************/
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- FaceList[nv] is a linked list of all Face structures for   -*/
  /*- which the smallest of the three vertex indices is nv.      -*/
  /*--------------------------------------------------------------*/
  SWGFace **FaceLists=(SWGFace **)mallocEC(NumVertices*sizeof(SWGFace *));
  memset(FaceLists,0,NumVertices*sizeof(SWGFace *));

  /*--------------------------------------------------------------*/
  /*- VertexUsed[nv] = true if vertex # nv is a vertex of any face*/
  /*- in the geometry. (Used below in the determination of the    */
  /*- number of "interior" vertices.)                             */
  /*--------------------------------------------------------------*/
  bool *VertexUsed=(bool *)mallocEC(NumVertices*sizeof(bool));
  memset(VertexUsed,0,NumVertices*sizeof(bool));

  /*--------------------------------------------------------------*/
  /*- EVNumEdges[nv] is the number of exterior edges connected to */
  /*- vertex nv (only filled in for nv=exterior vertex). This     */
  /*- number should be 2 for all exterior vertices.               */
  /*- EVEdges[nv][0] and EVEdges[nv][1] are pointers to the two   */
  /*- RWGEdge structures connected to vertex nv (again only       */
  /*- used for nv=exterior vertex)                                */
  /*--------------------------------------------------------------*/
#if 0
  EVNumEdges=(int *)mallocEC(NumVertices*sizeof(RWGEdge));
  memset(EVNumEdges,0,NumVertices*sizeof(int));
  EVEdges=(RWGEdge ***)mallocEC(NumVertices*sizeof(RWGEdge **));
  EVEdges[0]=(RWGEdge **)mallocEC(2*NumVertices*sizeof(RWGEdge *)); 
  for(nv=1; nv<NumVertices; nv++)
   EVEdges[nv]=EVEdges[nv-1]+2;
#endif

  /****************************************************************/
  /* construct a list of SWGFace structures by going over every   */
  /* face of every tetrahedron:                                   */ 
  /*                                                              */
  /*  a. if the face does not exist in the list (we have not      */
  /*     visited this edge before) then add it to the list        */
  /*  b. if the face does exist in the list, then we are visiting */
  /*     it for the second time, so it is an interior face, and   */
  /*     we mark it as such.                                      */
  /*  c. if we find ourselves visiting a face more than two       */
  /*     times then the topology of the mesh is defective.        */
  /*                                                              */
  /* note that when we first create a new SWGFace structure for a */
  /* face, we assign the Tet it came from as its 'positive'       */
  /* Tet, but we don't yet assign it an 'Index' in the overall    */
  /* problem, because it may turn out to be an exterior face. we  */
  /* also set the iQM, iMTet, and MIndex fields to -1 in view     */
  /* of this possibility.                                         */
  /*                                                              */
  /* if we subsequently come across the face for the second time  */
  /* (because it turned out to be attached to a second tet)       */
  /* we call this second tet the 'negative' tet associated to     */
  /* the face, and now we do assign an Index to the face, as well */
  /* as fill in its iQM, iMFace, and MIndex fields.               */
  /*                                                              */
  /* thus, at the conclusion of this loop, any Face structure that*/
  /* still has Index==-1 is an exterior face. such a structure    */
  /* will have valid data stored for its iV1, iV2, iV3, iQP,      */
  /* PTet, iPTet, and PIndex fields, but not its iQM, MTet, iMTet,*/
  /* or MIndex fields.                                            */
  /****************************************************************/
  NumInteriorFaces=NumTotalFaces=0;
  for(int nt=0; nt<NumTets; nt++)
   for(int nf=0; nf<4; nf++)   /* loop over tetrahedron faces */
    { 
      SWGTet *T=Tets[nt];

      /***************************************************************/
      /* get and sort the indices of the three vertices of face #nf. */
      /* Note: face #nf is the face opposite vertex #nf.             */
      /***************************************************************/
      int VI[3];
      VI[0]  = T->VI[ (nf+1)%4 ];
      VI[1]  = T->VI[ (nf+2)%4 ];
      VI[2]  = T->VI[ (nf+3)%4 ];

      int iMin = 0;
      if ( VI[1] < VI[iMin] ) iMin=1;
      if ( VI[2] < VI[iMin] ) iMin=2;
      int iMax = 0;
      if ( VI[1] > VI[iMax] ) iMax=1;
      if ( VI[2] > VI[iMax] ) iMax=2;

      int iVMin = VI[iMin];
      int iVMid = VI[3-iMin-iMax];
      int iVMax = VI[iMax];

      int iQ = T->VI[nf]; // vertex opposite face #nf

      VertexUsed[iVMin]=VertexUsed[iVMid]=VertexUsed[iVMax]=true;

      /**********************************************************************/
      /* look for this face in list of faces connected to vertex # iVMin    */
      /**********************************************************************/
      SWGFace *F;
      for(F=FaceLists[iVMin]; F; F=F->Next)
       if ( F->iV2==iVMid && F->iV3==iVMax )
        break;
      
      if ( F )
       { 
         /***************************************************************/
         /* we have encountered this face twice before ******************/
         /***************************************************************/
         if ( F->iMTet != -1 )
          ErrExit("%s: invalid mesh topology: face %i of tet %i also belongs to tets %i and %i ",
                      MeshFileName,nf,nt,F->iPTet,F->iMTet);

         /***************************************************************/
         /* we have encountered this face once before                   */
         /***************************************************************/
         F->iQM    = iQ;
         F->iMTet  = T->Index;
         F->MIndex = nf;
         F->Index  = NumInteriorFaces++;

         /* bounding radius is max distance from centroid to any vertex */
         F->Radius=VecDistance(F->Centroid, Vertices + 3*(F->iQP) );
         F->Radius=fmax(F->Radius, VecDistance(F->Centroid,Vertices + 3*F->iV1));
         F->Radius=fmax(F->Radius, VecDistance(F->Centroid,Vertices + 3*F->iV2));
         F->Radius=fmax(F->Radius, VecDistance(F->Centroid,Vertices + 3*F->iV3));
         F->Radius=fmax(F->Radius, VecDistance(F->Centroid,Vertices + 3*F->iQM));
       }
      else
       { 
         /******************************************************************/
         /* we are encountering this face for the first time. create a     */ 
         /* new SWGFace structure for this face and chain it in to the     */
         /* linked list of SWGFace structures connected to vertex iVMin.   */
         /******************************************************************/
         NumTotalFaces++;
         F=(SWGFace *)mallocEC(sizeof *F);
         F->Next=FaceLists[iVMin];
         FaceLists[iVMin]=F;

         F->iV1=iVMin;
         F->iV2=iVMid;
         F->iV3=iVMax;
         F->iQP=iQ;
         F->iQM=-1;

         F->iPTet  = T->Index;
         F->iMTet  = -1;

         F->PIndex = nf;
         F->MIndex = -1;

         F->Index = -1;

         InitSWGFace(F);
       };
    };

  /*--------------------------------------------------------------*/
  /*- now go back through our list of all faces:                  */
  /*-  a. put each SWGFace structure corresponding to an interior */
  /*      face or exterior face into the Faces array. Note that   */
  /*-     the first NumInteriorFaces slots in this array          */
  /*      correspond to interior faces, and all remaining slots   */
  /*      correspond to exterior faces.                           */
  /*-  c. While we're at it, update the 'FI' ('face index') fields*/
  /*-     of all Tet structures to indicate the indices that we   */
  /*-     have now assigned the various faces.                    */
  /*--------------------------------------------------------------*/
  Faces=(SWGFace **)mallocEC(NumTotalFaces*sizeof(Faces[0]));
  int NumExteriorFaces=0;
  for(int nv=0; nv<NumVertices; nv++)
   for(SWGFace *F=FaceLists[nv]; F; F=F->Next)
    { 
      if (F->Index==-1)  /* F is an exterior face */
       { 
         F->Index=NumInteriorFaces + NumExteriorFaces;
         Faces[ F->Index ] = F;
         NumExteriorFaces++;
       }
      else               /* F is an interior face */
       Faces[F->Index]=F;

      Tets[ F->iPTet ]->FI[ F->PIndex ] = F->Index;
      if ( F->iMTet !=-1 )
       Tets[ F->iMTet ]->FI[ F->MIndex ] = F->Index;
   };
  if ( (NumInteriorFaces + NumExteriorFaces) != NumTotalFaces )
   ErrExit("%s:%i: internal error (%i!=%i)",__FILE__,__LINE__,
            NumInteriorFaces+NumExteriorFaces, NumTotalFaces);

  /*--------------------------------------------------------------*/
  /*- count unused vertices --------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumUnusedVertices=0;
  for(int nv=0; nv<NumVertices; nv++)
   if (!VertexUsed[nv]) 
    NumUnusedVertices++;

  /*--------------------------------------------------------------*/
  /*- deallocate temporary storage -------------------------------*/
  /*--------------------------------------------------------------*/
  free(FaceLists);
  free(VertexUsed);

}

} // namespace buff
