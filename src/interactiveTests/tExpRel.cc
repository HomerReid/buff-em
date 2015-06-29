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
 * tExpRel.cc
 *
 * homer reid       -- 6/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fenv.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libbuff.h"

using namespace scuff;
using namespace buff;


#define II cdouble(0.0,1.0)

namespace buff{

void ExpRel23(cdouble x, cdouble *ExpRel2, cdouble *ExpRel3);

              }

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  double Theta=0.0;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Theta",              PA_DOUBLE,  1, 1, (void *)&Theta,          0, "theta angle"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  Theta *= M_PI / 180.0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Zeta = exp(II*Theta);
  FILE *f=fopen("tExpRel.out","w");
  for(double zMag = 1.0e-12; zMag<=100.0; zMag*=exp(0.1*log(10.0)))
   {
     cdouble z = zMag * Zeta;

     cdouble ER2Series, ER3Series;

     ExpRel23(z, &ER2Series, &ER3Series);

     cdouble ER2BF = exp(z) - 1.0 - z;
     cdouble ER3BF = ER2BF - 0.5*z*z;

     fprintf(f,"%e %e ",zMag,Theta);
     fprintf(f,"%e %e ",real(ER2Series),imag(ER2Series));
     fprintf(f,"%e %e ",real(ER3Series),imag(ER3Series));
     fprintf(f,"%e %e ",real(ER2BF),imag(ER2BF));
     fprintf(f,"%e %e ",real(ER3BF),imag(ER3BF));
     fprintf(f,"\n");
     
   };
  fclose(f);

}
