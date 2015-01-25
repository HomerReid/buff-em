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
 * ReadGMSHFile.cc -- subroutine of the SWGVolume class constructor
 *
 * homer reid    -- 3/2007 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;

namespace buff {

/*************************************************************/
/* constants needed in this file only  ***********************/
/*************************************************************/
#define TYPE_TETRAHEDRON 4
#define TYPE_POINT       15

#define NODE_START_KEYWORD1 "$NOD"
#define NODE_START_KEYWORD2 "$Nodes"

#define NODE_END_KEYWORD1   "$ENDNOD"
#define NODE_END_KEYWORD2   "$EndNodes"

#define ELM_START_KEYWORD1  "$ELM"
#define ELM_START_KEYWORD2  "$Elements"

#define ELM_END_KEYWORD1    "$ENDELM"
#define ELM_END_KEYWORD2    "$EndElements"

#define FORMAT_LEGACY 0
#define FORMAT_NEW    1

/*************************************************************/
/* Read vertices and panels from a GMSH .msh file to specify */
/* a surface.                                                */
/* If MeshTag==-1, then all panels are read.                 */
/* Otherwise, only panels on the specified physical region   */
/* are read.                                                 */
/*************************************************************/
void SWGVolume::ReadGMSHFile(FILE *MeshFile, const GTransformation *OTGT)
{
 
  /*------------------------------------------------------------*/
  /* Read lines until we hit the keyword indicating the start   */
  /* of the 'nodes' section, then read the number of nodes.     */
  /*------------------------------------------------------------*/
  char buffer[100];
  int LineNum=0;
  bool KeywordFound=false;
  char *MFN = MeshFileName;
  int WhichMeshFormat;
  while(!KeywordFound)
   { 
     if (!fgets(buffer,100,MeshFile))
      ErrExit("%s: failed to find node start keyword",MFN);
     LineNum++;
     if( !strncmp(buffer,NODE_START_KEYWORD1,strlen(NODE_START_KEYWORD1)))
      { WhichMeshFormat=FORMAT_LEGACY; KeywordFound=true; }
     else if( !strncmp(buffer,NODE_START_KEYWORD2,strlen(NODE_START_KEYWORD2)))
      { WhichMeshFormat=FORMAT_NEW; KeywordFound=true; }
   };
  if ( !fgets(buffer,100,MeshFile) || !(NumVertices=atoi(buffer)) )
   ErrExit("%s: invalid number of nodes",MFN);
  LineNum++;

  /*------------------------------------------------------------*/
  /*- Read in the vertices (which GMSH calls 'nodes.')          */
  /*- Note that the numbering of the vertices in GMSH does not  */
  /*- necessarily correspond to their ordering in the mesh      */
  /*- file. To remedy this situation, we construct a mapping    */
  /*- between GMSH's vertices indices and our internal vertex   */
  /*- indices, which works like this: The vertex that GMSH      */
  /*- calls 'node 3' is stored in slot GMSH2HR[3] within our    */
  /*- internal Vertices array.                                  */ 
  /*------------------------------------------------------------*/
  int NodeIndex;
  int *GMSH2HR=(int *)mallocEC(2*NumVertices*sizeof(int));
  for(int i=0; i<2*NumVertices; i++)
   GMSH2HR[i]=-1;
  Vertices=(double *)mallocEC(3*NumVertices*sizeof(double));
  for (int nv=0; nv<NumVertices; nv++)
   { if (!fgets(buffer,100,MeshFile))
      ErrExit("%s: too few nodes",MFN);
     LineNum++;
     int nConv=sscanf(buffer,"%i %le %le %le",&NodeIndex,
                          Vertices+3*nv,Vertices+3*nv+1,Vertices+3*nv+2); 
     if(nConv!=4)
      ErrExit("%s:%i: invalid node specification",MFN,LineNum); 
     if (NodeIndex>2*NumVertices)
      ErrExit("%s: internal error in ReadGMSHFile",MFN);
     GMSH2HR[NodeIndex]=nv;
   }; /* for (nv=0; nv<NumVertices; nv++) */
   
  /*------------------------------------------------------------*/
  /*- Apply one-time geometrical transformation (if any) to all */ 
  /*- vertices.                                                 */ 
  /*------------------------------------------------------------*/
  if (OTGT) OTGT->Apply(Vertices, NumVertices);
 
  /*------------------------------------------------------------*/
  /*- Eliminate any redundant vertices from the vertex list.   -*/
  /*- Explain me in greater detail please.                     -*/
  /*------------------------------------------------------------*/
  int NumRedundantVertices=0;
  for(int i=0; i<NumVertices; i++)
   for(int j=i+1; j<NumVertices; j++)
    if ( VecDistance(Vertices+3*i, Vertices+3*j) < 1.0e-6 )
     {
       /* remap all references to my node j so that they now refer to my node i*/
       for(int jGMSH=0; jGMSH<2*NumVertices; jGMSH++)
        if (GMSH2HR[jGMSH]==j)
         GMSH2HR[jGMSH]=i;

       NumRedundantVertices++;
     };
 
  /*------------------------------------------------------------*/
  /* Confirm that the next two lines in the file are the        */
  /* end-of-nodes-section keyword and the start-of-elements-section */
  /* keyword, then read the number of elements.                 */
  /*------------------------------------------------------------*/
  if ( !fgets(buffer,100,MeshFile) )
   ErrExit("%s: bad file format (nodes section not terminated)",MFN);
  LineNum++;

  if ( WhichMeshFormat==FORMAT_LEGACY ) 
   { if ( strncmp(buffer,NODE_END_KEYWORD1,strlen(NODE_END_KEYWORD1)))
      ErrExit("%s:%i: unexpected keyword",MFN,LineNum);
   }
  else
   { if ( strncmp(buffer,NODE_END_KEYWORD2,strlen(NODE_END_KEYWORD2)))
      ErrExit("%s:%i: unexpected keyword",MFN,LineNum);
   };

  if ( !fgets(buffer,100,MeshFile) )
   ErrExit("%s: bad file format (elements section not initiated)",MFN);
  LineNum++;

  if ( WhichMeshFormat==FORMAT_LEGACY ) 
   { if ( strncmp(buffer,ELM_START_KEYWORD1,strlen(ELM_START_KEYWORD1)))
      ErrExit("%s:%i: unexpected keyword",MFN,LineNum);
   }
  else
   { if ( strncmp(buffer,ELM_START_KEYWORD2,strlen(ELM_START_KEYWORD2)))
      ErrExit("%s:%i: unexpected keyword",MFN,LineNum);
   }

  if ( !fgets(buffer,100,MeshFile) )
   ErrExit("%s: bad file format (invalid number of elements)",MFN);
  LineNum++;

  int NumElements;
  int nConv=sscanf(buffer,"%i",&NumElements);
  if (nConv!=1 || NumElements<0) 
   ErrExit("%s:%i: invalid number of elements",MFN,LineNum);
 
  /*------------------------------------------------------------*/
  /*- Now read each line of the elements section.               */ 
  /*------------------------------------------------------------*/
  NumTets=0;
  Tets=(SWGTet **)mallocEC(NumElements * sizeof(Tets[0]));
  int VI[4];
  for (int ne=0; ne<NumElements; ne++)
   { 
     if (!fgets(buffer,100,MeshFile))
      ErrExit("too few elements in input file");
     LineNum++;

     int ElType;
     if (WhichMeshFormat==FORMAT_LEGACY)
      { 
        int ElNum, RegPhys, RegElem, NodeCnt;
        nConv=sscanf(buffer,"%i %i %i %i %i %i %i %i",
                     &ElNum,&ElType,&RegPhys,&RegElem,&NodeCnt,VI,VI+1,VI+2);
        if (nConv<5)
         ErrExit("%s:%i: invalid element specification",MFN,LineNum);
      }
     else
      { 
        int ElNum, nTags, nRead;
        nConv=sscanf(buffer,"%i %i %i%n",&ElNum,&ElType,&nTags,&nRead);
        int bufPos=nRead;
        if (nConv<3)
         ErrExit("%s:%i: invalid element specification",MFN,LineNum);

        // read the first 'tag,' which should be the physical region
        int RegPhys;
        if (nTags==0)
         RegPhys=0;
        else 
         sscanf(buffer+bufPos,"%i%n",&RegPhys,&nRead);
        bufPos+=nRead;

        // read any remaining tags 
        for(int nt=0; nt<nTags-1; nt++)
         { int iDummy;
           sscanf(buffer+bufPos,"%i%n",&iDummy,&nRead);
           bufPos+=nRead;
         };

        // finally, read the vertex indices. 
        nConv=sscanf(buffer+bufPos,"%i %i %i %i",VI,VI+1,VI+2,VI+3);

        if (ElType==TYPE_TETRAHEDRON && nConv!=4)
         ErrExit("%s:%i: invalid element specification",MFN,LineNum);
        else if (ElType==TYPE_POINT && nConv!=1)
         ErrExit("%s:%i: invalid element specification",MFN,LineNum);
      };
   
     /* we only process elements that are tetrahedra */
     switch(ElType)
      { 
        case TYPE_TETRAHEDRON:
          Tets[NumTets]=NewSWGTet(Vertices,
                                  GMSH2HR[ VI[0] ],
                                  GMSH2HR[ VI[1] ],
                                  GMSH2HR[ VI[2] ],
                                  GMSH2HR[ VI[3] ]
                                    );
          Tets[NumTets]->Index=NumTets;
          NumTets++;
          break;

        default:
          //ErrExit("%s:%i: unknown element type %i",MFN,LineNum,ElType);
          //fprintf(stderr,"%s:%i: warning: ignoring unknown element type %i\n",MF,LineNum,ElType);
          break;

      }; // switch(ElType)
   }; //for (ne=0; ne<NumElements; ne++)

  free(GMSH2HR);
  fclose(MeshFile);
 
} 

} // namespace buff
