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
 * TTaylorDuffy.h:  header file for the tetrahedron Taylor-Duffy method
 *                  for evaluating tetrahedra-product integrals over pairs
 *                  of tetrahedra with common vertices
 *
 * homer reid        1/2015
 */
#ifndef TTAYLORDUFFY_H
#define TTAYLORDUFFY_H

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <complex>

#include <libhrutil.h>
#include <libscuff.h>

using namespace scuff;

namespace buff { 

//values for the WhichK parameter
#define TTD_RP                   0
#define TTD_HELMHOLTZ            1
#define TTD_GRADHELMHOLTZ        2
#define NUMKS                    3

//values for the WhichP parameter
#define TTD_UNITY                0
#define TTD_BDOTBP               1

#define TTD_RXUNITY              2
#define TTD_RYUNITY              3
#define TTD_RZUNITY              4
#define TTD_TXUNITY              5
#define TTD_TYUNITY              6
#define TTD_TZUNITY              7

#define TTD_RXBDOTBP             8
#define TTD_RYBDOTBP             9
#define TTD_RZBDOTBP            10
#define TTD_TXBDOTBP            11
#define TTD_TYBDOTBP            12
#define TTD_TZBDOTBP            13

#define NUMPS                   14

//values for the WhichCase parameter (=number of common vertices)
#define TTD_COMMONTET            4
#define TTD_COMMONTRIANGLE       3
#define TTD_COMMONEDGE           2
#define TTD_COMMONVERTEX         1

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct TTDArgStruct
 { 
    // mandatory input fields 
    int WhichCase;
 
    int NumPKs;
    int *PIndex;
    int *KIndex;
    cdouble *KParam;

    double *V1, *V2, *V3, *V4;
    double *V2P, *V3P, *V4P;

    // output fields 
    cdouble *Result, *Error;
    int nCalls;
 
    // optional input fields
    double *Q, *QP;
    double XTorque[3];
    int RegionOnly;
    double AbsTol, RelTol;
    int MaxEval;

 } TTDArgStruct;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void TTaylorDuffy(TTDArgStruct *Args);
void InitTTDArgs(TTDArgStruct *Args);

} // namespace buff { 

#endif // TTAYLORDUFFY_H
