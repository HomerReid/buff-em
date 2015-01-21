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
 * buff-test-LFField.cc  -- buff-em unit test involving a dielectric 
 *                       -- sphere in a low-frequency external E-field
 *
 * homer reid            -- 1/2015
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
void DofprintVec(FILE *f, void *v, bool Complex, int Length,
                 const char *fmt, bool CR)
{
  char format[100];
  if (fmt)
   snprintf(format,100,"%s ",fmt);

  if (Complex)
   { cdouble *cv=(cdouble *)v;
     for(int n=0; n<Length; n++)
      fprintf(f,format,real(cv[n]),imag(cv[n]));
   }
  else
   { double *rv=(double *)v;
     for(int n=0; n<Length; n++)
      fprintf(f,format,rv[n]);
   };

  if (CR) 
   fprintf(f,"\n");
}

void fprintVec(FILE *f, double *v, int Length=3, const char *format="%+.8e")
{  DofprintVec(f, (void *)v, false, Length, format, false); }

void fprintVecCR(FILE *f, double *v, int Length=3, const char *format="%+.8e")
{  DofprintVec(f, (void *)v, false, Length, format, true); }

void fprintVec(FILE *f, cdouble *v, int Length=3, const char *format="%+.8e %+.8e")
{  DofprintVec(f, (void *)v, true, Length, format, false); }

void fprintVecCR(FILE *f, cdouble *v, int Length=3, const char *format="%+.8e %+.8e")
{  DofprintVec(f, (void *)v, true, Length, format, true); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SetLogFileName("buff-test-LFField.log");
  Log("buff-test-LFField running on %s",GetHostName());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SWGGeometry *G = new SWGGeometry("E10Sphere_533.buffgeo");
  HMatrix *M     = G->AllocateVIEMatrix();
  HVector *J     = G->AllocateRHSVector();
  cdouble Omega  = 0.1;

  HMatrix *XMatrix=new HMatrix("EPFile.XAxis",LHM_TEXT,"-ncol 3");
  if (XMatrix->ErrMsg)
   ErrExit(XMatrix->ErrMsg);

  cdouble E0[3]  = {1.0, 0.0, 0.0};
  double nHat[3] = {0.0, 0.0, 1.0};
  PlaneWave *PW  = new PlaneWave(E0, nHat);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  G->AssembleVIEMatrix(Omega, M);
  G->AssembleRHSVector(Omega, PW, J);
  M->LUFactorize();
  M->LUSolve(J);
  HMatrix *F1Matrix = G->GetFields(0,  J, Omega, XMatrix);
  HMatrix *F2Matrix = G->GetFields(PW, 0, Omega, XMatrix);

  FILE *f=fopen("buff-test-LFField.dat","w");
  for(int nx=0; nx<XMatrix->NR; nx++)
   {  
     double X[3];
     cdouble F1[6], F2[6];
     XMatrix->GetEntriesD(nx,"0:2",X);
     F1Matrix->GetEntries(nx,"0:5",F1);
     F2Matrix->GetEntries(nx,"0:5",F2);
     fprintVec(f,X,3);
     fprintVec(f,F1,6);
     fprintVecCR(f,F2,6);
   };
  fclose(f);
 
}
