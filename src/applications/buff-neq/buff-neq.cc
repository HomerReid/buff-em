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
 * buff-neq   -- a standalone code within the buff-em suite
 *             -- for implementing the fluctuating-surface-current
 *             -- approach to nonequilibrium phenomena (more 
 *             -- specifically, for computing heat radiation, 
 *             -- heat transfer, and nonequilibrium casimir forces) 
 *
 * homer reid  -- 5/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "buff-neq.h"
#include <libhrutil.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

#define MAXTEMPS 10    // max number of objects for which temperatures may be set

#define MAXSTR   1000

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SpeedTest(char *Greeting)
{
  double Elapsed;
  int Dim=1500;

  char *str=getenv("DIMDIM");
  if (str) sscanf(str,"%i",&Dim);
  printf("\n** %s: Dim=%i: \n",Greeting,Dim);

  HMatrix *Z1=new HMatrix(Dim,Dim,LHM_COMPLEX);
  HMatrix *Z2=new HMatrix(Dim,Dim,LHM_COMPLEX);

  srand48(time(0));
  for(int nr=0; nr<Dim; nr++)
   for(int nc=0; nc<Dim; nc++)
    Z1->SetEntry(nr,nc, drand48()+II*drand48() );

  Z2->Zero();
  for(int nr=0; nr<Dim; nr++)
   Z2->SetEntry(nr,nr,1.0);

  printf("** LU Factorize: ");
  Tic(); Z1->LUFactorize(); Elapsed=Toc(); 
  printf(" %.3g s\n",Elapsed);

  printf("** LU Invert (2): ");
  Tic(); Z1->LUSolve(Z2); Elapsed=Toc(); 
  printf(" %.3g s\n",Elapsed);

  printf("** LU Invert (1): ");
  Tic(); Z1->LUInvert(); Elapsed=Toc(); 
  printf(" %.3g s\n",Elapsed);

  double Norm=0.0, DNorm=0.0;
  cdouble z1, z2;

  for(int nr=0; nr<Dim; nr++)
   for(int nc=0; nc<Dim; nc++)
    { z1=Z1->GetEntry(nr,nc);
      z2=Z2->GetEntry(nr,nc);
      Norm+=abs(z1);
      DNorm+=abs(z1-z2);
    };
  printf("** Norm, DNorm, Ratio = %e, %e, %e\n",Norm,DNorm,DNorm/Norm);
  printf("\n\n");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *TransFile=0;

  int Power=0;
  int XForce=0;
  int YForce=0;
  int ZForce=0;
  int XTorque=0;
  int YTorque=0;
  int ZTorque=0;

  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
  double OmegaMin=0.00;              int nOmegaMin;
  double OmegaMax=-1.0;              int nOmegaMax;

  char *TempStrings[2*MAXTEMPS];     int nTempStrings;

  double AbsTol=0.0;
  double RelTol=5.0e-2;
  int Intervals=25;

  char *FileBase=0;

  bool UseExistingData=false;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,  0,             "list of geometrical transformation"},
/**/     
     {"Power",          PA_BOOL,    0, 1,       (void *)&Power,      0,             "compute power transfer"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,     0,             "compute X-force"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,     0,             "compute Y-force"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,     0,             "compute Z-force"},
     {"XTorque",        PA_BOOL,    0, 1,       (void *)&XTorque,    0,             "compute X-torque"},
     {"YTorque",        PA_BOOL,    0, 1,       (void *)&YTorque,    0,             "compute Y-torque"},
     {"ZTorque",        PA_BOOL,    0, 1,       (void *)&ZTorque,    0,             "compute Z-torque"},
/**/     
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,   "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles,  "list of (angular) frequencies"},
     {"OmegaMin",       PA_DOUBLE,  1, 1,       (void *)&OmegaMin,   &nOmegaMin,    "lower integration limit"},
     {"OmegaMax",       PA_DOUBLE,  1, 1,       (void *)&OmegaMax,   &nOmegaMax,    "upper integration limit"},
/**/     
     {"Temperature",    PA_STRING,  2, MAXTEMPS, (void *)TempStrings, &nTempStrings,  "set object xx temperature to xx"},
/**/     
     {"FileBase",       PA_STRING,  1, 1,       (void *)&FileBase,   0,             "prefix for names of .log, .flux, .out files"},
/**/     
     {"AbsTol",         PA_DOUBLE,  1, 1,       (void *)&AbsTol,     0,             "absolute tolerance for frequency quadrature"},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,     0,             "relative tolerance for frequency quadrature"},
     {"Intervals",      PA_INT,     1, 1,       (void *)&Intervals,  0,             "number of intervals for frequency quadrature"},
/**/
     {"UseExistingData", PA_BOOL,   0, 1,       (void *)&UseExistingData, 0,        "read existing data from .flux files"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  if (FileBase)
   SetLogFileName("%s.log",FileBase);
  else
   SetLogFileName("buff-neq.log");

  Log("buff-neq running on %s",GetHostName());

  /*******************************************************************/
  /* determine which output quantities were requested ****************/
  /*******************************************************************/
  int QuantityFlags=0;
  if (Power)  QuantityFlags|=QFLAG_POWER;
  if (XForce) QuantityFlags|=QFLAG_XFORCE;
  if (YForce) QuantityFlags|=QFLAG_YFORCE;
  if (ZForce) QuantityFlags|=QFLAG_ZFORCE;
  if (XTorque) QuantityFlags|=QFLAG_XTORQUE;
  if (YTorque) QuantityFlags|=QFLAG_YTORQUE;
  if (ZTorque) QuantityFlags|=QFLAG_ZTORQUE;
  if (QuantityFlags==0)
   ErrExit("you must specify at least one quantity to compute");

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run simulations                         */
  /*******************************************************************/
  HVector *OmegaPoints=0, *OmegaPoints0;
  int nFreq, nOV, NumFreqs=0;
  if (nOmegaFiles==1) // first process --OmegaFile option if present
   { 
     OmegaPoints=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaPoints->ErrMsg)
      ErrExit(OmegaPoints->ErrMsg);
     NumFreqs=OmegaPoints->N;
     Log("Read %i frequencies from file %s.",NumFreqs,OmegaFile);
   };

  // now add any individually specified --Omega options
  if (nOmegaVals>0)
   { 
     NumFreqs += nOmegaVals;
     OmegaPoints0=OmegaPoints;
     OmegaPoints=new HVector(NumFreqs, LHM_COMPLEX);
     nFreq=0;
     if (OmegaPoints0)
      { for(nFreq=0; nFreq<OmegaPoints0->N; nFreq++)
         OmegaPoints->SetEntry(nFreq, OmegaPoints0->GetEntry(nFreq));
        delete OmegaPoints0;
      };
     for(nOV=0; nOV<nOmegaVals; nOV++)
      OmegaPoints->SetEntry(nFreq+nOV, OmegaVals[nOV]);
     Log("Read %i frequencies from command line.",nOmegaVals);
   };

  /*******************************************************************/
  /* check that the user didn't simultaneously ask for a discret set */
  /* of frequencies and a frequency range over which to integrate;   */
  /* if a range was specified check that it makes sense              */
  /*******************************************************************/
  if ( OmegaPoints ) 
   { if ( nOmegaMin>0 || nOmegaMax>0 )
      ErrExit("--OmegaMin/--OmegaMax options may not be used with --Omega/--OmegaFile");
     if ( nTempStrings>0 )
      ErrExit("--Temperature option may not be used with --Omega/--OmegaFile");
     Log("Computing spectral density at %i frequencies.",NumFreqs);
   }
  else
   { 
     if ( nOmegaMin==1 && OmegaMin<0.0 )
      ErrExit("invalid value specified for --OmegaMin");
     if ( nOmegaMax==1 && OmegaMax<OmegaMin )
      ErrExit("invalid value specified for --OmegaMax");

     if ( OmegaMax==-1.0 )
      Log("Integrating over range Omega=(%g,infinity).",OmegaMin);
     else
      Log("Integrating over range Omega=(%g,%g).",OmegaMin,OmegaMax);
   };

  /*******************************************************************/
  /* create the BNEQData structure that contains all the info needed */
  /* to evaluate the neq transfer at a single frequency              */
  /*******************************************************************/
  BNEQData *BNEQD=CreateBNEQData(GeoFile, TransFile, QuantityFlags, FileBase);
  SWGGeometry *G=BNEQD->G;
  BNEQD->UseExistingData   = UseExistingData;

  /*******************************************************************/
  /* process any temperature specifications **************************/
  /*******************************************************************/
  double TEnvironment=0.0;
  double *TObjects=(double *)malloc(G->NumObjects*sizeof(double));
  memset(TObjects, 0, G->NumObjects*sizeof(double));
  if (nTempStrings)
   { 
     for(int nts=0; nts<nTempStrings; nts++)
      { 
        // first try to parse the number specified for the temperature 
        double TTemp;
        if ( 1!=sscanf(TempStrings[2*nts+1],"%le",&TTemp) )
         ErrExit("invalid temperature (%s) passed for --temperature option",TempStrings[2*nts+1]);

        // now try to figure out which temperature we need to set 
        if (!strcasecmp(TempStrings[2*nts],"ENVIRONMENT"))
         { TEnvironment=TTemp;
           Log("Setting environment temperature to %g kelvin.\n",TTemp);
           printf("Setting environment temperature to %g kelvin.\n",TTemp);
         }
        else 
         { int WhichObject;
           G->GetObjectByLabel(TempStrings[2*nts],&WhichObject);
           if (WhichObject==-1)
            ErrExit("unknown object (%s) passed for --temperature option",TempStrings[2*nts]);
           TObjects[WhichObject]=TTemp;
           Log("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
           printf("Setting temperature of object %s to %g kelvin.\n",TempStrings[2*nts],TTemp);
         };
      };
   };
         
  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  int OutputVectorLength = BNEQD->NumTransformations * G->NumObjects* G->NumObjects* BNEQD->NQ;
  double *I = new double[ OutputVectorLength ];
  if (NumFreqs>0)
   { for (nFreq=0; nFreq<NumFreqs; nFreq++)
      GetFlux(BNEQD, OmegaPoints->GetEntry(nFreq), I);
   }
  else
   { 
      double *E = new double[ OutputVectorLength ];
      EvaluateFrequencyIntegral2(BNEQD, OmegaMin, OmegaMax,
                                 TObjects, TEnvironment,
                                 Intervals, I, E);
      delete[] E;
   };
  delete[] I;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
