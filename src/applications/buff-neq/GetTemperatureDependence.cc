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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <libhrutil.h>
#include <libMDInterp.h>
#include <libhmat.h>
#include <libSGJC.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NUMPFT       8
//#define OMEGA_COLUMN 1
//#define DATA_COLUMN_OFFSET 3
#define OMEGA_COLUMN 0
#define DATA_COLUMN_OFFSET 2

#define HBAROMEGA02 9.491145534e-06
#define BOLTZMANNK 4.36763e-4
double GetThetaFactor(double Omega, double T)
{ 
  if (T==0.0)
   return 0.0;
  return Omega / ( exp( Omega/(BOLTZMANNK*T) ) - 1.0 );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct IntegrandData 
 {
   Interp1D *I1D;
   double Temperature;
   double TEnvironment;
   FILE *LogFile;
 } IntegrandData;

int MyIntegrand(unsigned ndim, const double *x, void *UserData,
                unsigned fdim, double *fval)
{ 
  (void) ndim;

  IntegrandData *Data = (IntegrandData *)UserData;
  Interp1D *I1D       = Data->I1D;
  double T            = Data->Temperature;
  double TEnvironment = Data->TEnvironment;
  FILE *LogFile       = Data->LogFile;

  double Omega=x[0];
  I1D->Evaluate(Omega, fval);
    
  double Factor = (   GetThetaFactor(Omega, T)
                    - GetThetaFactor(Omega, TEnvironment)
                  );
  Factor*=HBAROMEGA02;

  if (LogFile)
   { fprintf(LogFile,"%e %e ",Omega,Factor);
     for(unsigned nf=0; nf<fdim; nf++)
      fprintf(LogFile,"%e ",fval[nf]);
     fprintf(LogFile,"\n");
   };

  for(unsigned nf=0; nf<fdim; nf++)
   fval[nf] *= Factor;

  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *SIFluxFile=0;
  double TFile=0.0;
  double TMin=1.0;
  double TMax=300.0;
  double TStep=1.0;
  double TEnv=0.0;

  char *FileBase=0;

  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"SIFluxFile",     PA_STRING,  1, 1,       (void *)&SIFluxFile, 0, ".SIFlux file"},
/**/     
     {"FileTemperature", PA_DOUBLE,  1, 1,      (void *)&TFile,      0, "temperature at which .SIFlux file was written"},
/**/     
     {"TMin",           PA_DOUBLE,  1, 1,      (void *)&TMin,        0, "minimum output temperature"},
     {"TMax",           PA_DOUBLE,  1, 1,      (void *)&TMax,        0, "maximum output temperature"},
     {"TStep",          PA_DOUBLE,  1, 1,      (void *)&TStep,       0, "output temperature step"},
/**/
     {"TEnvironment",   PA_DOUBLE,  1, 1,      (void *)&TEnv,        0, "environment temperature"},
/**/
     {"FileBase",       PA_STRING,  1, 1,      (void *)&FileBase,    0, "base filename for output files"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (SIFluxFile==0)
   OSUsage(argv[0], OSArray, "--SIFluxFile option is mandatory");
 
  if (!FileBase)
   FileBase=strdup(GetFileBase(SIFluxFile));

  /***************************************************************/
  /* read the data file and sort in ascending order of frequency */
  /***************************************************************/
  HMatrix *SIFluxMatrix=new HMatrix(SIFluxFile);
  SIFluxMatrix->Sort(OMEGA_COLUMN);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NOmega=SIFluxMatrix->NR;
  double *OmegaValues=SIFluxMatrix->GetEntriesD(":", OMEGA_COLUMN);
  double *DataValues =(double *)mallocEC(NOmega*NUMPFT*sizeof(double));
  for(int nOmega=0; nOmega<NOmega; nOmega++)
   { double Omega = SIFluxMatrix->GetEntryD(nOmega,OMEGA_COLUMN);
     double Theta = (TFile==0.0) ? 1.0 : GetThetaFactor(Omega,TFile);
     double PreFactor = (TFile==0.0) ? 1.0 : HBAROMEGA02*Theta;
     for (int nq=0; nq<NUMPFT; nq++)
      { int nc=DATA_COLUMN_OFFSET+nq;
        double Data=SIFluxMatrix->GetEntryD(nOmega,nc);
        DataValues[nOmega*NUMPFT + nq]=Data/PreFactor;
      };
   };
  Interp1D *I1D=new Interp1D(OmegaValues, DataValues, NOmega, NUMPFT);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double OmegaMin = OmegaValues[0];
  double OmegaMax = OmegaValues[NOmega-1];
  struct IntegrandData ID, *Data=&ID;
  Data->I1D=I1D;
  Data->TEnvironment=TEnv;
  FILE *f=vfopen("%s.ByT","w",FileBase);
  for(double T=TMin; T<=TMax+0.5*TStep; T+=TStep)
   { 
     Data->Temperature=T;
 
     if (T==TMin)
      Data->LogFile=vfopen("%s.T%g.Log","w",FileBase,T);
     else if ( fabs(T-TMax) < 0.5*TStep )
      Data->LogFile=vfopen("%s.T%g.Log","w",FileBase,T);
     else
      Data->LogFile=0;

     int MaxEvals=0;
     double AbsTol=0.0; 
     double RelTol=1.0e-3;
     double Result[NUMPFT], Error[NUMPFT];
     pcubature(NUMPFT, MyIntegrand, (void *)Data, 1,
               &OmegaMin, &OmegaMax, MaxEvals, AbsTol, RelTol,
               ERROR_INDIVIDUAL, Result, Error);

     if (Data->LogFile)
      fclose(Data->LogFile);

     fprintf(f,"%e ",T);
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",Result[nq]);
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",Error[nq]);
     fprintf(f,"\n");
     fflush(f);
   };
  fclose(f);

  printf("Temperature data written to file %s.ByT.\n", GetFileBase(SIFluxFile));
  printf("Thank you for your support.\n");

}
