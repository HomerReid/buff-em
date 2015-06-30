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

void Mem(const char *s=0)
{
  if (s==0) s="default:";
  printf("\n*\n* %s\n*\n",s);
  unsigned long MU[7];
  GetMemoryUsage(MU);
  printf("%10lu %10lu %10lu %10lu\n",MU[0],MU[1],MU[2],MU[3]);
  printf("%10lu %10lu %10lu      \n",MU[4],MU[5],MU[6]);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetLogFileName("buff-test-FIBBICache.log");
  Log("buff-test-LFField running on %s",GetHostName());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry("E10Sphere_533.buffgeo");
  SWGVolume *O   = G->Objects[0];
  int nfA        = 0;
  int nfB        = 0;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double Data1[6], Data2[6];

  unsigned long M0=GetMemoryUsage();

  FIBBICache *GCache=new FIBBICache();
  printf("Alloc GCache: %lu\n",GetMemoryUsage()-M0);

  GCache->GetFIBBIData(O, nfA, O, nfB, Data1);
  printf("GetData 1   : %lu\n",GetMemoryUsage()-M0);

  GCache->GetFIBBIData(O, nfA, O, nfB, Data2);
  printf("GetData 2   : %lu\n",GetMemoryUsage()-M0);

  delete GCache;
  printf("After free  : %lu\n",GetMemoryUsage()-M0);

}
