/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * SVTensor.cc -- implementation of SVTensor class describing
 *             -- spatially-varying tensors
 *
 * homer reid  -- 12/2009
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <libhrutil.h>

#include "SVTensor.h"
#include "cmatheval.h"

#define MAXSTR 200

/***************************************************************/
/* constructor *************************************************/
/***************************************************************/
SVTensor::SVTensor(const char *FileName, bool IsMatProp)
{
   Name=strdupEC(FileName);
   ErrMsg=0;
   NumConstants=0;
   for(int nx=0; nx<3; nx++)
    for(int ny=0; ny<3; ny++)
     QExpression[nx][ny]=0;
   NumMPs=0;
   Homogeneous=false;
   Isotropic=false;

   /*--------------------------------------------------------------*/
   /*- detect cases in which the material is homogeneous and       */
   /*- described by a MatProp                                      */
   /*--------------------------------------------------------------*/
   if (IsMatProp)
    { 
      MPs[0]=new MatProp(Name);
      if ( MPs[0]->ErrMsg)
       ErrExit(MPs[0]->ErrMsg);
      
      Log("Created isotropic material with MatProp=%s",Name);
      MatProp::SetLengthUnit(1.0e-6);
      Homogeneous=Isotropic=true;
      NumMPs=1;
      return;
    };
 
   /*--------------------------------------------------------------*/
   /*- try to open the file ---------------------------------------*/
   /*--------------------------------------------------------------*/
   char *WhichDir;
   FILE *f=fopenPath(getenv("BUFF_SVTENSOR_PATH"),Name,"r",&WhichDir);
   if (!f)
    { ErrMsg = vstrdup("could not find file %s",Name);
      return;
    };
   Log("Found SVTensor file %s/%s",WhichDir,Name);

   /*--------------------------------------------------------------*/
   /*- parse the file ---------------------------------------------*/
   /*--------------------------------------------------------------*/
   ErrMsg = Parse(f);
   fclose(f); 
   if (ErrMsg) return;

   /*--------------------------------------------------------------*/
   /*- initialize cevaluators -------------------------------------*/
   /*--------------------------------------------------------------*/
   for(int nx=0; nx<3; nx++)
    for(int ny=0; ny<3; ny++)
     { 
       void *Expr = QExpression[nx][ny];
       if (Expr==0) continue;

       cevaluator_set_var_index(Expr, "w",     0);
       cevaluator_set_var_index(Expr, "x",     1);
       cevaluator_set_var_index(Expr, "y",     2);
       cevaluator_set_var_index(Expr, "z",     3);
       cevaluator_set_var_index(Expr, "r",     4);
       cevaluator_set_var_index(Expr, "Theta", 5);
       cevaluator_set_var_index(Expr, "Phi",   6);

       // add variables named Eps1, Eps2, ... for scuff-em matprops
       for(int nmp=0; nmp<NumMPs; nmp++)
        { char str[10];
          snprintf(str,10,"Eps%i",nmp+1);
          cevaluator_set_var_index(Expr, str, 7+nmp);
        };

       for(int nc=0; nc<NumConstants; nc++)
        cevaluator_set_var(Expr,ConstantNames[nc],ConstantValues[nc]);
     };
}

/*--------------------------------------------------------------*/
/*- if the function specification contains any strings of the   */
/*- form MAT_XX, then replace them with Eps1, Eps2, etc. and    */
/*- add scuff-em material properties for material XX            */
/*--------------------------------------------------------------*/
void SVTensor::ExtractMPs(char *str, char *FileName, int LineNum)
{ 
  while( char *p=strstr(str, "MP_") )
   {
     // get the portion of the string between _ and the next space
     // into MPName
     char MPName[MAXSTR], *q=p+3;
     int n=0;
     while( *q && *q!=' ' && n<(MAXSTR-1))
      { MPName[n++]=*q;
        *q++ = ' ';
      };
     MPName[n]=0;

     // see if we already have an MP with this name
     int MPIndex=0;
     for(int nmp=0; MPIndex==0 && nmp<NumMPs; nmp++)
      if ( !strcasecmp(MPs[nmp]->Name,MPName) )
       MPIndex=nmp+1;

     // if not, try to create one 
     if (MPIndex==0)
      { if (NumMPs==MAXMPS)
         ErrExit("%s:%i: ",FileName,LineNum,"too many MP_ statements");
        MPs[NumMPs] = new MatProp(MPName);
        if (MPs[NumMPs]->ErrMsg)
         ErrExit("%s:%i: ",FileName,LineNum,MPs[NumMPs]->ErrMsg);
        MPIndex=++NumMPs;
        MatProp::SetLengthUnit(1.0e-6);
      };

     // replace the string e.g. "MP_Gold" with e.g. "Eps1"
     snprintf(MPName,MAXSTR,"Eps%i",MPIndex);
     strncpy(p, MPName, strlen(MPName));

   };
 
}

/***************************************************************/
/* constructor helper function *********************************/
/***************************************************************/
char *SVTensor::Parse(FILE *f)
{
  /*--------------------------------------------------------------*/
  /*- parse lines one at a time ----------------------------------*/
  /*--------------------------------------------------------------*/
  int LineNum=0;
  char Line[MAXSTR];
  char *FileName=Name;
  while ( fgets(Line,MAXSTR,f) )
   { 
     LineNum++;

     /*--------------------------------------------------------------*/
     /*- skip blank lines and comments  -----------------------------*/
     /*--------------------------------------------------------------*/
     char *p=Line;
     while ( *p && isspace(*p) )
      p++;
     if ( *p==0 || *p=='#' )
      continue;

     /*--------------------------------------------------------------*/
     /*- strip off trailing semicolon and/or carriage return       --*/
     /*--------------------------------------------------------------*/
     int Len=strlen(p);
     char *pp;
     if ( (pp=strchr(p,';')) )
      *pp=0;
     if ( (pp=strchr(p,'\n')) )
      *pp=0;

     /*--------------------------------------------------------------*/
     /*- hand off to the MatProp constructor to parse any            */
     /*- MATERIAL...ENDMATERIAL sections                             */
     /*--------------------------------------------------------------*/
     if ( !strncmp(p,"MATERIAL",8) )
      {
        char *Tokens[2];
        int nTokens=Tokenize(p, Tokens, 2);

        if ( nTokens==1 )
         ErrExit("%s:%i: no name given for MATERIAL ",FileName,LineNum);
        else if ( nTokens>2 )
         ErrExit("%s:%i: syntax error",FileName,LineNum);
         
        ErrMsg=AddMaterialToMatPropDataBase(f, FileName, Tokens[1], &LineNum);
        if (ErrMsg)
         ErrExit("%s:%i: %s",FileName,LineNum,ErrMsg);

        continue;
      };

     /*--------------------------------------------------------------*/
     /*- detect lines of the form -----------------------------------*/
     /*-   Q   =                  -----------------------------------*/
     /*- Eps   =                  -----------------------------------*/
     /*- Qxx   =                  -----------------------------------*/
     /*- EpsXX =                  -----------------------------------*/
     /*--------------------------------------------------------------*/
     char *q=0;
     if ( Len>=3 && !strncasecmp(p,"Q",1) )
      q=p+1;
     else if ( Len>=5 && !strncasecmp(p,"Eps",1) )
      q=p+3;
     if ( q )
      { 
        pp=strchr(q,'=');
        int MN[2]={0,0};

        // isotropic case
        if ( *q=='=' || *q==' ' )
         Isotropic=true;
        // anisotropic case
        else 
         for(int n=0; n<2; n++)
          { // the user can specify either e.g. EpsXY or Eps12
            if( toupper(q[n])>='X' && toupper(q[n])<='Z')
             MN[n] = toupper(q[n]) - 'X';
            else if (q[n]>='1' && q[n]<='3')
             MN[n] = q[n] - '1';
            else 
             return vstrdup("%s:%i: invalid expression",FileName,LineNum);
          };

        int M=MN[0], N=MN[1];
        ExtractMPs(pp+1,FileName,LineNum);
        QExpression[M][N]=cevaluator_create(pp+1);
        if (!QExpression[M][N])
         return vstrdup("%s:%i: invalid expression",FileName,LineNum);
      }
     /*--------------------------------------------------------------*/
     /*- detect constant specifications of the form NAME = [number] -*/
     /*--------------------------------------------------------------*/
     else if ( pp )
      { 
        double CV;
        if ( 1!=sscanf(pp+1,"%le",&CV) )
         return vstrdup("%s:%i: invalid constant value",FileName,LineNum);

        if ( NumConstants==MAXCONSTANTS )
         return vstrdup("%s:%i: too many constants (max allowed is %i)",FileName,LineNum,MAXCONSTANTS);

        while ( isspace (*(pp-1) ) ) 
          pp--;
        *pp=0; 
        ConstantNames[NumConstants]=strdup(p);
        ConstantValues[NumConstants]=CV;
        NumConstants++;
      }
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     else
      return vstrdup("%s:%i: syntax error",FileName,LineNum);

    }; // while ( fgets(Line,MAXSTR,f) )

  return 0;

}

/***************************************************************/
/* destructor **************************************************/
/***************************************************************/
SVTensor::~SVTensor()
{
  free(Name);

  for(int nx=0; nx<3; nx++)
   for(int ny=0; ny<3; ny++)
    if ( QExpression[nx][ny] )
     cevaluator_destroy(QExpression[nx][ny]);

  for(int nmp=0; nmp<NumMPs; nmp++)
   delete MPs[nmp];

  for(int nc=0; nc<NumConstants; nc++)
   free(ConstantNames[nc]);

}  

/***************************************************************/
/* get tensor components at a given frequency and location  ****/
/***************************************************************/
void SVTensor::Evaluate(cdouble Omega, double x[3], cdouble Q[3][3])
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Homogeneous && Isotropic)
   { cdouble QMP = MPs[0]->GetEps(Omega);
     for(int nx=0; nx<3; nx++)
      for(int ny=0; ny<3; ny++)
       Q[nx][ny] = (nx==ny) ? QMP : 0.0;
     return;
   };

  /***************************************************************/
  /* elsewhere I wrote things so that the maximum number of MP_  */
  /* inclusions is the variable MAXMPS, but here I have hard-    */
  /* coded it to 3                                               */
  /***************************************************************/
  if (MAXMPS!=3) ErrExit("%s:%i: internal error",__FILE__,__LINE__);
  int NumVars=10;
  static const char *VNames[10] = {"w",
                                   "x", "y", "z",
                                   "r", "Theta", "Phi",
                                   "Eps1", "Eps2", "Eps3"};
  cdouble VValues[10];

  VValues[0] = Omega;
  VValues[1] = x[0];
  VValues[2] = x[1];
  VValues[3] = x[2];
  VValues[4] = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
  VValues[5] = atan2( sqrt(x[0]*x[0] + x[1]*x[1]), x[2] );
  VValues[6] = atan2( x[1], x[2] );
  for(int nmp=0; nmp<NumMPs; nmp++)
   VValues[7 + nmp] = MPs[nmp]->GetEps(Omega);
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0; nx<3; nx++)
   for(int ny=0; ny<3; ny++)
    Q[nx][ny] = (nx==ny) ? 1.0 : 0.0; 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char **VVNames = const_cast<char **>(VNames);
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    if (QExpression[Mu][Nu])
     Q[Mu][Nu] = cevaluator_evaluate(QExpression[Mu][Nu], NumVars,
                                     VVNames, VValues);

  // enforce symmetry of off-diagonals
  for(int Mu=0; Mu<3; Mu++)
   { int Nu = (Mu+1)%3;
     cdouble QAvg = 0.5*(Q[Mu][Nu] + Q[Nu][Mu]);
     Q[Mu][Nu]=Q[Nu][Mu]=QAvg;
   };

  // if Qyy and Qzz weren't specified but Qxx was,
  // set Qyy=Qzz=Qxx
  if (QExpression[0][0])
   { if (QExpression[1][1]==0) Q[1][1]=Q[0][0];
     if (QExpression[2][2]==0) Q[2][2]=Q[0][0];
   };

}

cdouble SVTensor::Evaluate(cdouble Omega, double x[3])
{ cdouble Q[3][3];
  Evaluate(Omega, x, Q);
  return Q[0][0];
}

double SVTensor::EvaluateD(cdouble Omega, double x[3])
 { return real(Evaluate(Omega, x)); }
