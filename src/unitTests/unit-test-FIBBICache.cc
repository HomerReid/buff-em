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
 * buff-test-FIBBICache.cc  -- buff-em unit test for FIBBI cache
 *                          -- functionality
 *
 * homer reid               -- 6/2015
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libbuff.h"

using namespace scuff;
using namespace buff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  int NumTests=0, NumFailed=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetLogFileName("buff-test-FIBBICache.log");
  Log("buff-test-FIBBICache running on %s",GetHostName());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry("E10P1ISphere_48.buffgeo");
  SWGVolume *O   = G->Objects[0];
  int nfA        = O->NumInteriorFaces - 1;
  int nfB        = nfA;
  
  /***************************************************************/
  /* first test whether a single cache entry can be properly     */
  /* stored and retrieved                                        */
  /***************************************************************/
  unsigned long M0=GetMemoryUsage();

  FIBBICache *GCache=new FIBBICache();

  Tic(true);
  double Data1[6], Data2[6];
  GCache->GetFIBBIData(O, nfA, O, nfB, Data1);
  unsigned long M1;
  double Time1=Toc(&M1);
  Log("GetData 1 : %.3e us, %8lu b allocated", Time1, M1);

  Tic(true);
  GCache->GetFIBBIData(O, nfA, O, nfB, Data2);
  unsigned long M2;
  double Time2=Toc(&M2);
  Log("GetData 2 : %.3e us, %8lu b allocated", Time2, M2);
  
  NumTests++;
  if ( memcmp(Data2, Data1, 6*sizeof(double) ) )
   { Log(" Data 2 != Data 1");
     NumFailed++;
   };

  delete GCache;
  Log("After cache delete: %lu",GetMemoryUsage()-M0);

  /***************************************************************/
  /* next test whether the cache file is properly written to     */
  /* disk by the VIE matrix assembly                             */
  /***************************************************************/
  unlink("Sphere_48.GCache");
  HMatrix *M=G->AllocateVIEMatrix();
  Log("Assembling VIE matrix the first time...");
  G->AssembleVIEMatrix(1.0, M);
  Log("After VIE matrix assembly: %8lu MB allocated",
       GetMemoryUsage() / (1<<20));

  Tic(true);
  double Data3[6];
  GCache=G->ObjectGCaches[0];
  int Size1=GCache->Size();
  GCache->GetFIBBIData(O, nfA, O, nfB, Data3);
  unsigned long M3;
  double Time3=Toc(&M3);
  Log("GetData 3 : %.3e us, %8lu b allocated", Time3, M3);

  NumTests++;
  if ( memcmp(Data3, Data1, 6*sizeof(double) ) )
   { Log(" Data 3 != Data 1");
     NumFailed++;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Log("Reading cache from disk...");
  GCache=new FIBBICache(G->Objects[0]->MeshFileName);
  int Size2=GCache->Size();
  NumTests++;
  if (Size2!=Size1)
   { Log(" read only %i records (should have been %i)!",Size2,Size1);
     NumFailed++;
   };

  Tic(true);
  double Data4[6];
  GCache->GetFIBBIData(O, nfA, O, nfB, Data4);
  unsigned long M4;
  double Time4=Toc(&M4);
  Log("GetData 4 : %.3e us, %8lu b allocated", Time4, M4);
  NumTests++;
  if ( memcmp(Data4, Data1, 6*sizeof(double) ) )
   { Log(" Data 4 != Data 1");
     NumFailed++;
   };

  Log("%i/%i tests passed.",NumTests-NumFailed,NumTests);
  return NumFailed;
}
