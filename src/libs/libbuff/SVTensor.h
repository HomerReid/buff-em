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
 * SVTensor.h  -- a simple C++ class for managing
 *             -- Spatially Varying Tensors
 *
 * homer reid  -- 5/2014
 */

#ifndef SVTENSOR_H
#define SVTENSOR_H

#define MAXMPS 3
#define MAXCONSTANTS 25

#include <libMatProp.h>

/***************************************************************/
/* MatProp class definition ************************************/
/***************************************************************/
class SVTensor
 { 
  public:  

   /***************************************************************/
   /* public class methods ****************************************/
   /***************************************************************/
   /* constructor */
   SVTensor(const char *FileName, bool IsMatProp=false);

  /* destructor */
   ~SVTensor();

   /* get the tensor components at a given frequency and position */
   void Evaluate(cdouble Omega, double x[3], cdouble Q[3][3]);

   /* if ErrMsg is not NULL after the class constructor is invoked, there  */
   /* was an error.                                                        */
   char *ErrMsg;

 // private:

   /***************************************************************/
   /***************************************************************/
   /***************************************************************/
   cdouble EvalExpression(void *Expression, cdouble Omega, double x[3]);

   //constructor helper functions
   void ExtractMPs(char *s, char *FileName, int LineNum);
   char *Parse(FILE *f);

   /***************************************************************/
   /* class data **************************************************/
   /***************************************************************/
   char *Name;
   MatProp *MPs[MAXMPS];      // pointers to SCUFF-EM mat props
   int NumMPs;
   bool Homogeneous;
   bool Isotropic;

   // opaque pointers to cmatheval parsed expressions for
   // components of Q
   void *QExpression[3][3];

   int NumConstants;
   char *ConstantNames[MAXCONSTANTS];
   double ConstantValues[MAXCONSTANTS];

 };

#endif
