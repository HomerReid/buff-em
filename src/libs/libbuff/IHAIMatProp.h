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
 * IHAIMatProp.h  -- a simple C++ class for managing material
 *                -- properties for InHomogeneous, AnIsotropic
 *                -- materials
 *
 * homer reid    -- 5/2014
 */

#ifndef IHAIMATPROP_H
#define IHAIMATPROP_H

#define MAXCONSTANTS 25

/***************************************************************/
/* MatProp class definition ************************************/
/***************************************************************/
class IHAIMatProp
 { 
  public:  

   /***************************************************************/
   /* public class methods ****************************************/
   /***************************************************************/
   /* constructor */
   IHAIMatProp(const char *IHAIMatFileName);

  /* destructor */
   ~IHAIMatProp();

   /* get epsilon and mu at a given frequency and position */
   void GetEpsMu(cdouble Omega, double x[3], 
                 cdouble Eps[3][3], cdouble Mu[3][3]);
   void GetEps(cdouble Omega, double x[3], cdouble Eps[3][3]);

   /* if ErrMsg is not NULL after the class constructor is invoked, there  */
   /* was an error.                                                        */
   char *ErrMsg;

 // private:

   /***************************************************************/
   /***************************************************************/
   /***************************************************************/
   cdouble EvalExpression(void *Expression, cdouble Omega, double x[3]);
   char *Parse(FILE *f);

   /***************************************************************/
   /* class data **************************************************/
   /***************************************************************/
   char *Name;
   cdouble ConstEps; // =0.0 for parsed materials 

   // opaque pointers to cmatheval parsed expressions for
   // components of epsilon and mu
   void *EpsExpression[3][3], *MuExpression[3][3];

   int NumConstants;
   char *ConstantNames[MAXCONSTANTS];
   double ConstantValues[MAXCONSTANTS];

 };

#endif
