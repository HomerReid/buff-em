/* Copyright (C) 2005-2014 M. T. Homer Reid
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
 * buff-cas3D.cc  -- computing the Casimir energy, force, and/or
 *                -- torque for a collection of interacting objects
 *
 * homer reid     -- 5/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "buff-cas3D.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXFREQ  10    // max number of frequencies 
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

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
//
  bool Energy=false;
  bool XForce=false;
  bool YForce=false;
  bool ZForce=false;
  double TorqueAxes[9];                     int nTorque;
  bool AllTorque=false;
//
  double Temperature=0.0;                   int nTemperature;
  double XiVals[MAXFREQ];	            int nXiVals;
  char *XiFile=0;
//
  int MaxXiPoints=10000;
  double AbsTol=0.0;
  double RelTol=1.0e-2;
  int Intervals=50;
//
  char *OutputFile=0;
  char *ByXiFile=0;
  char *LogFile=0;
//
  char *Cache=0;
  char *ReadCache[MAXCACHE];                int nReadCache;
  char *WriteCache=0;
//
  bool UseExistingData=false;
//
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"Geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,       0,             "geometry file"},
     {"TransFile",      PA_STRING,  1, 1,       (void *)&TransFile,     0,             "list of geometrical transformation"},
//
     {"Energy",         PA_BOOL,    0, 1,       (void *)&Energy,        0,             "compute Casimir energy"},
     {"XForce",         PA_BOOL,    0, 1,       (void *)&XForce,        0,             "compute x-directed Casimir force"},
     {"YForce",         PA_BOOL,    0, 1,       (void *)&YForce,        0,             "compute y-directed Casimir force"},
     {"ZForce",         PA_BOOL,    0, 1,       (void *)&ZForce,        0,             "compute z-directed Casimir force"},
     {"Torque",         PA_DOUBLE,  3, 3,       (void *)TorqueAxes,     &nTorque,      "compute Casimir torque about a given axis"},
     {"AllTorque",      PA_BOOL,    0, 1,       (void *)&AllTorque,     0,             "compute all three Casimir torque components"},
//
     {"Temperature",    PA_DOUBLE,  1, 1,       (void *)&Temperature,   &nTemperature, "temperature in Kelvin"},
     {"Xi",             PA_DOUBLE,  1, MAXFREQ, (void *)XiVals,         &nXiVals,      "imaginary frequency"},
     {"XiFile",         PA_STRING,  1, 1,       (void *)&XiFile,        0,             "list of --Xi values"},
//
     {"OutputFile",     PA_STRING,  1, 1,       (void *)&OutputFile,    0,             "name of frequency-integrated output file"},
     {"ByXiFile",       PA_STRING,  1, 1,       (void *)&ByXiFile,      0,             "name of frequency-resolved output file"},
     {"LogFile",        PA_STRING,  1, 1,       (void *)&LogFile,       0,             "name of log file"},
//
     {"MaxXiPoints",    PA_INT,     1, 1,       (void *)&MaxXiPoints,   0,             "maximum number of Xi integrand evaluations "},
     {"AbsTol",         PA_DOUBLE,  1, 1,       (void *)&AbsTol,        0,             "absolute tolerance for sums and integrations"},
     {"RelTol",         PA_DOUBLE,  1, 1,       (void *)&RelTol,        0,             "relative tolerance for sums and integrations"},
     {"Intervals",      PA_INT,     1, 1,       (void *)&Intervals,     0,             "number of subintervals for frequency quadrature"},
//
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,         0,             "read/write cache"},
     {"ReadCache",      PA_STRING,  1, MAXCACHE,(void *)ReadCache,      &nReadCache,   "read cache"},
     {"WriteCache",     PA_STRING,  1, 1,       (void *)&WriteCache,    0,             "write cache"},
//
     {"UseExistingData", PA_BOOL,   0, 1,       (void *)&UseExistingData, 0,           "use alternative method for energy calculation"},
//
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (LogFile)
   SetLogFileName(LogFile);
  else
   SetLogFileName("buff-cas3D.log");
  Log("buff-cas3D running on %s",GetHostName());

  /***************************************************************/
  /* try to create the geometry  *********************************/
  /***************************************************************/
  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");
  SWGGeometry *G = new SWGGeometry(GeoFile);

  /***************************************************************/
  /* process frequency-related options                           */
  /***************************************************************/
  HMatrix *XiPoints=0;
  if ( XiFile )
   { if (nXiVals>0)    ErrExit("--XiFile and --Xi options are mutually exclusive");
     if (nTemperature) ErrExit("--XiFile and --Temperature options are mutually exclusive");
     XiPoints = new HMatrix(XiFile,LHM_TEXT,"--nc 1 --strict");
     if (XiPoints->ErrMsg)
      ErrExit(XiPoints->ErrMsg);
     Log("Read %i Xi points from file %s.",XiPoints->NR, XiFile);
   }
  else if ( nXiVals>0 )
   { if (nTemperature) ErrExit("--Xi and --Temperature options are mutually exclusive");
     XiPoints = new HMatrix(nXiVals, 1);
     for(int nxv=0; nxv<nXiVals; nxv++) 
      XiPoints->SetEntry(nxv, 0, XiVals[nxv]);
     Log("Performing Casimir calculations at %i command-line Xi points.",nXiVals);
   }
  else if ( nTemperature==1 )
   Log("Computing full Matsubara-summed Casimir quantities at T=%e Kelvin.",Temperature);
  else
   Log("Computing full zero-temperature Casimir quantities.");

  /*******************************************************************/
  /* figure out which quantities to compute **************************/
  /*******************************************************************/
  int WhichQuantities=0, NumQuantities=0;
  if (Energy)    { NumQuantities++; WhichQuantities |= QUANTITY_ENERGY; };
  if (XForce)    { NumQuantities++; WhichQuantities |= QUANTITY_XFORCE; };
  if (YForce)    { NumQuantities++; WhichQuantities |= QUANTITY_YFORCE; };
  if (ZForce)    { NumQuantities++; WhichQuantities |= QUANTITY_ZFORCE; };
  if (nTorque>0) { NumQuantities++; WhichQuantities |= QUANTITY_TORQUE1; };
  if (nTorque>1) { NumQuantities++; WhichQuantities |= QUANTITY_TORQUE2; };
  if (nTorque>2) { NumQuantities++; WhichQuantities |= QUANTITY_TORQUE3; };
  if (AllTorque)
   { if (nTorque>0) ErrExit("--AllTorque and --Torque options are mutually exclusive");
     NumQuantities+=3; 
     WhichQuantities |= (QUANTITY_TORQUE1 + QUANTITY_TORQUE2 + QUANTITY_TORQUE3);
     nTorque=3;
     memset(TorqueAxes, 0, 9*sizeof(double));
     TorqueAxes[0]=TorqueAxes[4]=TorqueAxes[8]=1.0;
   };

  /*******************************************************************/
  /* create the BC3Data structure that contains all the info needed  */
  /* to evaluate the contributions of a single frequency and kBloch  */
  /* point to the Casimir quantities                                 */
  /*******************************************************************/
  BC3Data *BC3D=CreateBC3Data(G, TransFile, WhichQuantities, NumQuantities,
                              nTorque, TorqueAxes);

  BC3D->ByXiFileName       = ByXiFile; 
  BC3D->AbsTol             = AbsTol;
  BC3D->RelTol             = RelTol;
  BC3D->MaxXiPoints        = MaxXiPoints;
  BC3D->UseExistingData    = UseExistingData;

  if (BC3D->ByXiFileName==0)
   { BC3D->ByXiFileName=vstrdup("%s.byXi",GetFileBase(G->GeoFileName));
     FILE *f=fopen(BC3D->ByXiFileName,"a"); 
     if (!f) ErrExit("could not open file %s",BC3D->ByXiFileName);
     WriteFilePreamble(f, BC3D, PREAMBLE_BYXI);
     fclose(f);
   };

  /*******************************************************************/
  /* now switch off based on the requested frequency behavior to     */
  /* perform the actual calculations                                 */
  /*******************************************************************/
  double *EFT = new double[BC3D->NTNQ];
  double *Error=0; 
  if ( XiPoints )
   { 
     for (int nr=0; nr<XiPoints->NR; nr++)
      GetXiIntegrand(BC3D, XiPoints->GetEntryD(nr,0), EFT);
   }
  else if ( Temperature > 0.0)
   { 
     Error = new double[BC3D->NTNQ];
     GetMatsubaraSum(BC3D, Temperature, EFT, Error);
   }
  else
   { 
     Error = new double[BC3D->NTNQ];
     if (Intervals==0)
      GetXiIntegral_Adaptive(BC3D, EFT, Error);
     else
      GetXiIntegral_Fixed(BC3D, Intervals, EFT, Error);
   };

  /***************************************************************/
  /* write output file if we computed summed or integrated quantities */
  /***************************************************************/
  if (Error)
   { if (!OutputFile)
      OutputFile=vstrdup("%s.out",GetFileBase(G->GeoFileName));
     FILE *OutFile=CreateUniqueFile(OutputFile,1);
     WriteFilePreamble(OutFile, BC3D, PREAMBLE_OUT);
     int ntnq=0;
     for(int nt=0; nt<BC3D->NumTransformations; nt++)
      { fprintf(OutFile,"%s ",BC3D->GTCList[nt]->Tag);
        for(int nq=0; nq<BC3D->NumQuantities; nq++, ntnq++)
         fprintf(OutFile,"%e %e ",EFT[ntnq],Error[ntnq]);
        fprintf(OutFile,"\n");
      };
     fclose(OutFile); 
     delete[] Error;
   };

  delete[] EFT;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
