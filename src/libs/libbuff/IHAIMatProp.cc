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
 * IHAIMatProp.cc -- implementation of IHAIMatProp class
 *
 * homer reid     -- 12/2009
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <libhrutil.h>

#include "IHAIMatProp.h"
#include "cmatheval.h"

#define MAXSTR       1000

/***************************************************************/
/* constructor *************************************************/
/***************************************************************/
IHAIMatProp::IHAIMatProp(const char *IHAIMatFileName )
 {
   Name=strdupEC(IHAIMatFileName);
   ErrMsg=0;
   NumConstants=0;
   for(int nx=0; nx<3; nx++)
    for(int ny=0; ny<3; ny++)
     EpsExpression[nx][ny]=MuExpression[nx][ny]=0;
   ConstEps=0.0;
   

   /*--------------------------------------------------------------*/
   /*- detect filenames of the form CONST_EPS_xx-------------------*/
   /*--------------------------------------------------------------*/
   if (!strncasecmp(Name,"CONST_EPS_",10))
    { int Status=S2CD(Name+10,&ConstEps);
      if (Status!=0 || ConstEps==0.0)
       { ErrMsg = vstrdup("invalid material specification %s",Name);
         return;
       };
      Log("Created constant-epsilon material with Eps=%s",z2s(ConstEps));
      return;
    };
 
   /*--------------------------------------------------------------*/
   /*- try to open the file ---------------------------------------*/
   /*--------------------------------------------------------------*/
   FILE *f=fopen(Name,"r");
   if (!f)
    { ErrMsg = vstrdup("could not open file %s",IHAIMatFileName);
      return;
    };

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
     { void *Expr = EpsExpression[nx][ny];
       if (Expr==0) continue;
       cevaluator_set_var_index(Expr, "x", 0);
       cevaluator_set_var_index(Expr, "y", 1);
       cevaluator_set_var_index(Expr, "z", 2);
       cevaluator_set_var_index(Expr, "w", 3);
       for(int nc=0; nc<NumConstants; nc++)
        cevaluator_set_var(Expr,ConstantNames[nc],ConstantValues[nc]);
     };
}

/***************************************************************/
/* constructor helper function *********************************/
/***************************************************************/
char *IHAIMatProp::Parse(FILE *f)
{
  /*--------------------------------------------------------------*/
  /*- parse lines one at a time ----------------------------------*/
  /*--------------------------------------------------------------*/
  int LineNum=0;
  char Line[MAXSTR];
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

     int Len=strlen(p);
     char *pp=strchr(p,'=');
     if (pp) // strip off trailing semicolon and/or carriage return
      { char *ppp;
        if ( (ppp=strchr(pp,';')) )
         *ppp=0;
        if ( (ppp=strchr(pp,'\n')) )
         *ppp=0;
      };

     /*--------------------------------------------------------------*/
     /*- detect lines of the form -----------------------------------*/
     /*-  Eps(w,x,y,z) = ...      -----------------------------------*/
     /*--------------------------------------------------------------*/
     if ( pp && Len>=12 && !strncasecmp(p,"Eps(w,x,y,z)",12 ) )
      { EpsExpression[0][0] = cevaluator_create(pp+1);
        if (!EpsExpression[0][0])
         return vstrdup("%s:%i: invalid expression",Name,LineNum);
      }
     /*--------------------------------------------------------------*/
     /*- detect lines of the form -----------------------------------*/
     /*-  EpsMN(w,x,y,z) = ...    -----------------------------------*/
     /*--------------------------------------------------------------*/
     else if ( pp && Len>=14 && !strncasecmp(p,  "Eps",3 )
               && !strncasecmp(p+5,"(w,x,y,z)",9 )
             )
      { 
        // the '-1' here allows users to use one-based indexing
        // for tensor indices
        int M=p[3] - '0' - 1;
        int N=p[4] - '0' - 1;

        if ( M<0 || M>3 || N<0 || N>3)
         return vstrdup("%s:%i: invalid indices %c%c",Name,LineNum,p[3],p[4]);

        EpsExpression[M][N]=cevaluator_create(pp+1);
        if (!EpsExpression[M][N])
         return vstrdup("%s:%i: invalid expression",Name,LineNum);
      }
     /*--------------------------------------------------------------*/
     /*- detect constant specifications of the form NAME = [number] -*/
     /*--------------------------------------------------------------*/
     else if ( pp )
      { 
        double CV;
        if ( 1!=sscanf(pp+1,"%le",&CV) )
         return vstrdup("%s:%i: invalid constant value",Name,LineNum);

        if ( NumConstants==MAXCONSTANTS )
         return vstrdup("%s:%i: too many constants (max allowed is %i)",Name,LineNum,MAXCONSTANTS);

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
      return vstrdup("%s:%i: syntax error",Name,LineNum);
    }; // while ( fgets(Line,MAXSTR,f) )

  return 0;

}

/***************************************************************/
/* destructor **************************************************/
/***************************************************************/
IHAIMatProp::~IHAIMatProp()
{
  free(Name);

  for(int nx=0; nx<3; nx++)
   for(int ny=0; ny<3; ny++)
    { if ( EpsExpression[nx][ny] )
       cevaluator_destroy(EpsExpression[nx][ny]);
      if ( MuExpression[nx][ny] )
       cevaluator_destroy(MuExpression[nx][ny]);
    };

  for(int nc=0; nc<NumConstants; nc++)
   free(ConstantNames[nc]);

}  

/***************************************************************/
/* get eps and mu at a given frequency and location  ***********/
/***************************************************************/
void IHAIMatProp::GetEpsMu(cdouble Omega, double x[3],
                           cdouble Eps[3][3], cdouble Mu[3][3])
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (ConstEps!=0.0)
   { for(int nx=0; nx<3; nx++)
      for(int ny=0; ny<3; ny++)
       { Eps[nx][ny] = (nx==ny) ? ConstEps : 0.0;
          Mu[nx][ny] = (nx==ny) ? 1.0      : 0.0;
       };
     return;
   };

  static char *VNames[4] = {"x", "y", "z", "w"};
  cdouble VValues[4];

  VValues[0] = x[0];
  VValues[1] = x[1];
  VValues[2] = x[2];
  VValues[3] = Omega;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nx=0; nx<3; nx++)
   for(int ny=0; ny<3; ny++)
    Eps[nx][ny] = Mu[nx][ny] = (nx==ny) ? 1.0 : 0.0; 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( EpsExpression[0][0] )
   Eps[0][0]=Eps[1][1]=Eps[2][2]
    =cevaluator_evaluate(EpsExpression[0][0], 4, VNames, VValues);

  if ( EpsExpression[1][1] )
   Eps[1][1]
    =cevaluator_evaluate(EpsExpression[1][1], 4, VNames, VValues);

  if ( EpsExpression[2][2] )
   Eps[2][2]
    =cevaluator_evaluate(EpsExpression[2][2], 4, VNames, VValues);

  if ( EpsExpression[0][1] )
   { Eps[0][1]=cevaluator_evaluate(EpsExpression[0][1], 4, VNames, VValues);
     if (!EpsExpression[1][0]) Eps[1][0]=Eps[0][1];
   };

  if ( EpsExpression[1][0] )
   { Eps[1][0]=cevaluator_evaluate(EpsExpression[1][0], 4, VNames, VValues);
     if (!EpsExpression[0][1]) Eps[0][1]=Eps[1][0];
   };

  if ( EpsExpression[0][2] )
   { Eps[0][2]=cevaluator_evaluate(EpsExpression[0][2], 4, VNames, VValues);
     if (!EpsExpression[2][0]) Eps[2][0]=Eps[0][2];
   };

  if ( EpsExpression[2][0] )
   { Eps[2][0]=cevaluator_evaluate(EpsExpression[2][0], 4, VNames, VValues);
     if (!EpsExpression[0][2]) Eps[0][2]=Eps[2][0];
   };

  if ( EpsExpression[1][2] )
   { Eps[1][2]=cevaluator_evaluate(EpsExpression[1][2], 4, VNames, VValues);
     if (!EpsExpression[2][1]) Eps[2][1]=Eps[1][2];
   };

  if ( EpsExpression[2][1] )
   { Eps[2][1]=cevaluator_evaluate(EpsExpression[2][1], 4, VNames, VValues);
     if (!EpsExpression[1][2]) Eps[1][2]=Eps[2][1];
   };
}

void IHAIMatProp::GetEps(cdouble Omega, double x[3],
                         cdouble Eps[3][3])
{ cdouble Mu[3][3];
  GetEpsMu(Omega, x, Eps, Mu);
}
