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
 * tFaceInt
 *
 * homer reid       -- 6/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;
using namespace buff;

#define II cdouble(0.0,1.0)
#define NFUN 12

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData
 { double *Centroid;
 } IntegrandData;
 
void Integrand(double *x, double *b, double Divb,
               void *UserData, double *I)
{
  IntegrandData *Data = (IntegrandData *)UserData;
  double *Centroid = Data->Centroid;

  double R[3];
  R[0] = (x[0]-Centroid[0]);
  R[1] = (x[1]-Centroid[1]);
  R[2] = (x[2]-Centroid[2]);

  I[0]  = b[0];
  I[1]  = b[1];
  I[2]  = b[2];
  I[3]  = b[0]*R[0];
  I[4]  = b[0]*R[1];
  I[5]  = b[0]*R[2];
  I[6]  = b[1]*R[0];
  I[7]  = b[1]*R[1];
  I[8]  = b[1]*R[2];
  I[9]  = b[2]*R[0];
  I[10] = b[2]*R[1];
  I[11] = b[2]*R[2];
  
}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  int nf=-1;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING,  1, 1, (void *)&GeoFile,  0, ".buffgeo file"},
     {"nf",                 PA_INT,     1, 1, (void *)&nf,       0, "nf"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry(GeoFile);
  SWGVolume *O = G->Objects[0];
  
  srand48(time(0));
  if (nf==-1)
   nf = lrand48() % O->NumInteriorFaces;
  printf("--nf %i \n",nf);

  /***************************************************************/
  /* get moments using exact formulas  ***************************/
  /***************************************************************/
  double J[3], Q[3][3];
  GetDQMoments(O, nf, J, Q);

  /***************************************************************/
  /* get moments using numerical cubature ************************/
  /***************************************************************/
  IntegrandData MyData, *Data=&MyData;
  Data->Centroid = O->Faces[nf]->Centroid;
  double JQ[NFUN], JQError[NFUN];
  BFInt(O, nf, Integrand, (void *)Data,
        NFUN, JQ, JQError, 0, 0, 1.0e-8);

  /***************************************************************/
  /* compare *****************************************************/
  /***************************************************************/
  for(int Mu=0; Mu<3; Mu++)
   printf("J%i:  %+.8e   %+.8e   %.1e \n", 
           Mu, J[Mu],JQ[Mu],RD(J[Mu],JQ[Mu]));

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    printf("Q%i%i: %+.8e   %+.8e   %.1e (%e) \n", Mu, Nu, 
            Q[Mu][Nu],JQ[3 + 3*Mu + Nu],RD(Q[Mu][Nu],JQ[3 + 3*Mu + Nu]),
            Q[Mu][Nu]/JQ[3 + 3*Mu + Nu]
          );

}
