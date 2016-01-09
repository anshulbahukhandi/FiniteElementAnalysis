/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_CONSTANTS_H__
#define __RAJAN_CONSTANTS_H__

const int PDIM      = 2;        // problem dimensionality
const int NUMENODES = 4;        // number of element nodes
const int DOFPN     = 2;        // number of unknowns per node
const int NUMEP     = 4;        // number of element properties
const int NUMMATP   = 1;        // number of material properties
const int NUMSCOMP  = 3;        // number of element strains/stresses
                                // points per element
const double TOL = 1.0e-6;      // tolerance used in equation solver
const int MAXCHARS = 80;        // max. chars in problem description

#endif