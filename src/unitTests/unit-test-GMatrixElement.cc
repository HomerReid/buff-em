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
 * buff-test-TetTetInt.cc -- test of the three methods implemented in
 *                        -- buff-EM for computing 6-dimensional
 *                        -- integrals over tetrahedron-product domains
 *
 * homer reid             -- 1/2015
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libbuff.h"
#include "TetTetInt.h"

using namespace scuff;
using namespace buff;

#define FDIM 4

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(cdouble *V1, cdouble *V2, int N,
             const char *str1, const char *str2)
{ 
  //char FStr[10];
  //snprintf(FStr,10,"%+.%ie",Precision);

  printf(" n | %-25s | %-25s | RD      | Ratio\n",str1,str2);
  for(int n=0; n<N; n++)
   printf("%2i | (%+.4e,%+.4e) | (%+.4e,%+.4e) | %.1e | %.1e\n",n,
    real(V1[n]),imag(V1[n]), real(V2[n]),imag(V2[n]),
    RD(V1[n],V2[n]), abs(V1[n]/V2[n]));
  printf("\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(cdouble *V1, cdouble *V2, cdouble *V3, int N,
             const char *str1, const char *str2, const char *str3)
{ 
  //char FStr[10];
  //snprintf(FStr,10,"%+.%ie",Precision);

  printf(" n | %-25s | %-25s | RD      | %-25s | RD\n",str1,str2,str3);
  for(int n=0; n<N; n++)
   printf("%2i | (%+.4e,%+.4e) | (%+.4e,%+.4e) | %.1e | (%+.4e, %+.4e) | %.1e\n",n,
    real(V1[n]),imag(V1[n]),
    real(V2[n]),imag(V2[n]),
    RD(V1[n],V2[n]), 
    real(V3[n]),imag(V3[n]),
    RD(V1[n],V3[n]));
  printf("\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetLogFileName("buff-test-TetTetInt.log");
  Log("buff-test-TetTetInt running on %s",GetHostName());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry("TwoTetrahedra.buffgeo");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TetTetIntArgs MyTTIArgs, *TTIArgs=&MyTTIArgs;
  InitTetTetIntArgs(TTIArgs);
  TTIArgs->noa=0;
  TTIArgs->nta=0;
  TTIArgs->iQa=0;
  TTIArgs->nob=0;
  TTIArgs->ntb=0;
  TTIArgs->iQb=0;

  int WhichPoly[FDIM]  = { BUFF_POLY_UNITY,
                           BUFF_POLY_UNITY,
                           BUFF_POLY_UNITY,
                           BUFF_POLY_UNITY
                         };
  int WhichKernel[FDIM] = { BUFF_KERN_PSI,
                            BUFF_KERN_PSI,
                            BUFF_KERN_PSI,
                            BUFF_KERN_PSI,
                          };
  cdouble KParam[FDIM]  =  {1.0, 2.0, 3.0, 4.0};

  TTIArgs->fdim        = FDIM; 
  TTIArgs->WhichPoly   = WhichPoly;
  TTIArgs->WhichKernel = WhichKernel;
  TTIArgs->KParam      = KParam;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble IBF[FDIM];
  TTIArgs->CubatureMethod = BUFF_CUBATURE_BF;
  Tic();
  MyTetTetInt(G, TTIArgs, IBF);
  printf("BF time: %e ms\n",1.0e3*Toc());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble ITD[FDIM];
  TTIArgs->CubatureMethod = BUFF_CUBATURE_TD;
  Tic();
  MyTetTetInt(G, TTIArgs, ITD);
  printf("TD time: %e ms\n",1.0e3*Toc());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble ISI[FDIM];
  TTIArgs->CubatureMethod = BUFF_CUBATURE_SI;
  Tic();
  MyTetTetInt(G, TTIArgs, ISI);
  printf("SI time: %e ms\n",1.0e3*Toc());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Compare(ITD, ISI, IBF, FDIM, "TD", "SI", "BF");

}
