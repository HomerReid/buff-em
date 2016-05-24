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
 * buff-scatter.cc  -- a standalone code within the buff-EM suite for
 *                  -- solving problems involving the scattering of
 *                  -- electromagnetic radiation from an arbitrary
 *                  -- compact object
 *
 * homer reid       -- 6/2014
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "buff-scatter.h"

using namespace scuff;
using namespace buff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXPW    1     // max number of plane waves
#define MAXGB    1     // max number of gaussian beams
#define MAXPS    10    // max number of point sources
#define MAXFREQ  10    // max number of frequencies
#define MAXEPF   10    // max number of evaluation-point files
#define MAXCACHE 10    // max number of cache files for preload

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
  char *GeoFile=0;
//
  double pwDir[3*MAXPW];             int npwDir;
  cdouble pwPol[3*MAXPW];            int npwPol;
//
  double gbDir[3*MAXGB];             int ngbDir;
  cdouble gbPol[3*MAXGB];            int ngbPol;
  double gbCenter[3*MAXGB];          int ngbCenter;
  double gbWaist[MAXGB];             int ngbWaist;
//
  double psLoc[3*MAXPS];             int npsLoc;
  cdouble psStrength[3*MAXPS];       int npsStrength;
//
  char *IFFile=0;
//
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile;                   int nOmegaFiles;
//
  char *EPFiles[MAXEPF];             int nEPFiles;
//
  char *JPlotFile=0;
//
  char *PFTFile          = 0;
  char *OPFTFile         = 0;
  char *EMTPFTFile       = 0;
  char *MomentPFTFile    = 0;
  char *DSIPFTFile       = 0;
  double DSIRadius       = 10.0;
  int DSIPoints          = 302;
  int DSIPoints2         = 0;
  char *DSIMesh          = 0;
//
  char *MomentFile=0;

  int ExportMatrix=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,            "geometry file"},
/**/
     {"Omega",          PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,   &nOmegaVals,  "(angular) frequency"},
     {"OmegaFile",      PA_STRING,  1, 1,       (void *)&OmegaFile,  &nOmegaFiles, "list of (angular) frequencies"},
/**/
     {"pwDirection",    PA_DOUBLE,  3, MAXPW,   (void *)pwDir,       &npwDir,      "plane wave direction"},
     {"pwPolarization", PA_CDOUBLE, 3, MAXPW,   (void *)pwPol,       &npwPol,      "plane wave polarization"},
/**/
     {"gbDirection",    PA_DOUBLE,  3, MAXGB,   (void *)gbDir,       &ngbDir,      "gaussian beam direction"},
     {"gbPolarization", PA_CDOUBLE, 3, MAXGB,   (void *)gbPol,       &ngbPol,      "gaussian beam polarization"},
     {"gbCenter",       PA_DOUBLE,  3, MAXGB,   (void *)gbCenter,    &ngbCenter,   "gaussian beam center"},
     {"gbWaist",        PA_DOUBLE,  1, MAXGB,   (void *)gbWaist,     &ngbWaist,    "gaussian beam waist"},
/**/
     {"psLocation",     PA_DOUBLE,  3, MAXPS,   (void *)psLoc,       &npsLoc,      "point source location"},
     {"psStrength",     PA_CDOUBLE, 3, MAXPS,   (void *)psStrength,  &npsStrength, "point source strength"},
/**/
     {"IFFile",         PA_STRING,  1, 1,       (void *)&IFFile,     0,             "list of incident fields"},
/**/
     {"EPFile",         PA_STRING,  1, MAXEPF,  (void *)EPFiles,     &nEPFiles,    "list of evaluation points"},
/**/
     {"PFTFile",        PA_STRING,  1, 1,       (void *)&PFTFile,    0,            "name of PFT output file (computed by default EMTPFT method)"},
     {"EMTPFTFile",     PA_STRING,  1, 1,       (void *)&EMTPFTFile, 0,            "name of J \\dot E PFT output file"},
     {"OPFTFile",       PA_STRING,  1, 1,       (void *)&OPFTFile,   0,            "name of overlap PFT output file"},
     {"MomentPFTFile",  PA_STRING,  1, 1,       (void *)&MomentPFTFile,0,          "name of moment PFT output file"},
     {"DSIPFTFile",     PA_STRING,  1, 1,       (void *)&DSIPFTFile, 0,            "name of DSIPFT output file"},
     {"DSIMesh",        PA_STRING,  1, 1,       (void *)&DSIMesh,    0,            "mesh file for surface-integral PFT"},
     {"DSIRadius",      PA_DOUBLE,  1, 1,       (void *)&DSIRadius,  0,            "radius of bounding sphere for surface-integral PFT"},
     {"DSIPoints",      PA_INT,     1, 1,       (void *)&DSIPoints,  0,            "number of quadrature points for surface-integral PFT (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)"},
     {"DSIPoints2",     PA_INT,     1, 1,       (void *)&DSIPoints2, 0,            "number of quadrature points for DSIPFT second opinion"},
/**/
     {"JPlotFile",      PA_STRING,  1, 1,       (void *)&JPlotFile,  0,            "name of J-plot file"},
/**/
     {"MomentFile",     PA_STRING,  1, 1,       (void *)&MomentFile, 0,            "name of induced-dipole-moment output file"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaList=0;
  if (nOmegaFiles==1) // first process --OmegaFile option if present
   { 
     OmegaList=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaList->ErrMsg)
      ErrExit(OmegaList->ErrMsg);
   };

  // now add any individually specified --Omega options
  if (nOmegaVals>0)
   { 
     int N1 = (OmegaList ? OmegaList->N : 0);
     int N2 = nOmegaVals;
     HVector *NewOmegaList=new HVector(N1+N2, LHM_COMPLEX); 
     if (OmegaList)
      { memcpy(NewOmegaList->ZV, OmegaList->ZV, N1*sizeof(cdouble));
        delete OmegaList;
      };
     OmegaList=NewOmegaList;
     for(int n2=0; n2<N2; n2++)
      OmegaList->SetEntry(N1+n2, OmegaVals[n2]);
   };

  if ( !OmegaList || OmegaList->N==0)
   OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem   */
  /*******************************************************************/
  if ( npwPol != npwDir )
   ErrExit("numbers of --pwPolarization and --pwDirection options must agree");
  if ( ngbPol != ngbDir || ngbDir!=ngbCenter || ngbCenter!=ngbWaist )
   ErrExit("numbers of --gbPolarization, --gbDirection, --gbCenter, and --gbWaist options must agree ");
  if ( npsLoc!=npsStrength )
   ErrExit("numbers of --psLocation and --psStrength options must agree");

  IncField *IF=0;
  for(int npw=0; npw<npwPol; npw++)
   { PlaneWave *PW=new PlaneWave(pwPol + 3*npw, pwDir + 3*npw);
     PW->Next=IF;
     IF=PW;
   };
  for(int ngb=0; ngb<ngbCenter; ngb++)
   { GaussianBeam *GB=new GaussianBeam(gbCenter + 3*ngb, gbDir + 3*ngb, gbPol + 3*ngb, gbWaist[ngb]);
     GB->Next=IF;
     IF=GB;
   };
  for(int nps=0; nps<npsLoc; nps++)
   { PointSource *PS=new PointSource(psLoc + 3*nps, psStrength + 3*nps);
     PS->Next=IF;
     IF=PS;
   };

  IncFieldList *IFList=0;
  if (IFList!=0 && IF!=0)
   ErrExit("--IFFile is incompatible with other incident-field specifications");
  else if (IFFile!=0 && IF==0)
   IFList = ReadIncFieldList(IFFile);
  else if (IFFile==0 && IF!=0)
   IFList = AddIncFieldToList(IF,const_cast<char *>("Default"));

  /*******************************************************************/
  /* sanity check to make sure the user specified an incident field  */
  /* if one is required for the outputs the user requested           */
  /*******************************************************************/
  bool NeedIncidentField = (ExportMatrix==false);
  if ( NeedIncidentField && IFList==0 )
   ErrExit("you must specify at least one incident field source");

  /*******************************************************************/
  /* PFT options *****************************************************/
  /*******************************************************************/
  PFTOptions MyPFTOptions, *pftOptions = &MyPFTOptions;
  BUFF_InitPFTOptions(pftOptions);
  pftOptions->DSIMesh       = DSIMesh;
  pftOptions->DSIRadius     = DSIRadius;
  pftOptions->DSIPoints     = DSIPoints;

  char *DSIPFTFile2=0;
  if (DSIPFTFile && DSIPoints2!=0)
   DSIPFTFile2=vstrdup("%s.DSI%i",GetFileBase(DSIPFTFile),DSIPoints2);

  /*******************************************************************/
  /* create the BSData structure containing everything we need to    */
  /* execute scattering calculations                                 */
  /*******************************************************************/
  BSData MyBSData, *BSD=&MyBSData;

  SWGGeometry *G      = BSD->G   = new SWGGeometry(GeoFile);
  HMatrix *M          = BSD->M   = G->AllocateVIEMatrix();
                        BSD->RHS = G->AllocateRHSVector();
  HVector *J          = BSD->J   = G->AllocateRHSVector();
  BSD->IF             = 0;
  BSD->IFLabel        = 0;
  BSD->TransformLabel = 0;

  char GeoFileBase[MAXSTR];
  strncpy(GeoFileBase, GetFileBase(GeoFile), MAXSTR);

  /*******************************************************************/
  /* loop over frequencies *******************************************/
  /*******************************************************************/   
  char OmegaStr[MAXSTR];
  cdouble Omega;
  cdouble Eps, Mu;
  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
   { 
     Omega = OmegaList->GetEntry(nFreq);
     BSD->Omega = Omega;
     z2s(Omega, OmegaStr);
     Log("Working at frequency %s...",OmegaStr);

     /*******************************************************************/
     /* assemble VIE matrix at this frequency                           */
     /*******************************************************************/
     G->AssembleVIEMatrix(Omega, M);

     /*******************************************************************/
     /* export VIE matrix to a binary file if that was requested        */
     /*******************************************************************/
     if (ExportMatrix)
      { void *pCC=HMatrix::OpenMATLABContext("%s_%s",GeoFileBase,OmegaStr);
        M->ExportToMATLAB(pCC,"M");
        HMatrix::CloseMATLABContext(pCC);
      };

     /*******************************************************************/
     /* if the user requested no output options (for example, if she   **/
     /* just wanted to export the matrix to a binary file), don't      **/
     /* bother LU-factorizing the matrix or assembling the RHS vector. **/
     /*******************************************************************/
     if ( !NeedIncidentField )
      continue;

     /*******************************************************************/
     /* LU-factorize the VIE matrix to prepare for solving scattering   */
     /* problems                                                        */
     /*******************************************************************/
     Log("  LU-factorizing VIE matrix...");
     M->LUFactorize();

     /***************************************************************/
     /* loop over incident fields                                   */
     /***************************************************************/
     for(int nIF=0; nIF<IFList->NumIFs; nIF++)
      { 
        IF = BSD->IF = IFList->IFs[nIF];
        BSD->IFLabel = IFFile ? IFList->Labels[nIF] : 0;
        if (BSD->IFLabel)
         Log("  Processing incident field %s...",BSD->IFLabel);

        char IFStr[100]="";
        if (BSD->IFLabel)
         snprintf(IFStr,100,"_%s",BSD->IFLabel);

        /***************************************************************/
        /* set up the incident field profile and assemble the RHS vector */
        /***************************************************************/
        Log("  Assembling the RHS vector..."); 
        G->AssembleRHSVector(Omega, IF, J);
        if (JPlotFile)
         G->PlotCurrentDistribution(JPlotFile, J, "RHS_%s%s",CD2S(Omega),IFStr);
        BSD->RHS->Copy(J); // save a copy of the RHS vector for later

        /***************************************************************/
        /* solve the VIE system ****************************************/
        /***************************************************************/
        Log("  Solving the VIE system...");
        M->LUSolve(J);

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        if (JPlotFile)
         G->PlotCurrentDistribution(JPlotFile, J, "J_%s",CD2S(Omega),IFStr);

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        if (PFTFile)
         WritePFTFile(BSD, PFTFile, pftOptions, SCUFF_PFT_EMT);
   
        if (EMTPFTFile)
         WritePFTFile(BSD, EMTPFTFile, pftOptions, SCUFF_PFT_EMT);
   
        if (OPFTFile)
         WritePFTFile(BSD, OPFTFile, pftOptions, SCUFF_PFT_OVERLAP);
   
        if (MomentPFTFile)
         WritePFTFile(BSD, MomentPFTFile, pftOptions, SCUFF_PFT_MOMENTS);
   
        if (DSIPFTFile)
         { 
           pftOptions->DSIPoints = DSIPoints;
           WritePFTFile(BSD, DSIPFTFile, pftOptions, SCUFF_PFT_DSI);
   
           if (DSIPoints2)
            { pftOptions->DSIPoints = DSIPoints2;
              WritePFTFile(BSD, DSIPFTFile2, pftOptions, SCUFF_PFT_DSI);
            };
   
         };
   
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        if (MomentFile)
         WriteMomentFile(BSD, MomentFile);
   
        /*--------------------------------------------------------------*/
        /*- scattered fields at user-specified points ------------------*/
        /*--------------------------------------------------------------*/
        for(int nepf=0; nepf<nEPFiles; nepf++)
         ProcessEPFile(BSD, EPFiles[nepf]);

      }; // for(int nIF=0; nIF<IFList->NumIFs; nIF++)
   
   }; //  for(nFreq=0; nFreq<OmegaList->N; nFreqs++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}

