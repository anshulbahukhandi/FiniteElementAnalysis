/*********************************************
Utility Library Function
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_GETINTERACTIVE_H__
#define __RAJAN_GETINTERACTIVE_H__

#include <string>
#include <string>
#include <math.h>
#include <limits.h>
#include <float.h>

// int version
void GetInteractive (const std::string& szString, int& V);
void GetInteractive (const std::string& szString, int& V,
                     int nL, int nU);
void GetInteractive (const std::string& szString, int V[],
                     int nSize);
// long version
void GetInteractive (const std::string& szString, long& V);
void GetInteractive (const std::string& szString, long& V,
                     long lL, long lU);
// float version
void GetInteractive (const std::string& szString, float& V);
void GetInteractive (const std::string& szString, float& V,
                     float fL, float fU);
void GetInteractive (const std::string& szString, float V[],
                     int nSize);
// double version
void GetInteractive (const std::string& szString, double& V);
void GetInteractive (const std::string& szString, double& V,
                     double dL, double dU);
// string version
void GetInteractive (const std::string& szString, std::string& V,
                     int n);

// raw access and conversion functions
int GetIntValue (const std::string& szInput, int& nV);
int GetLongValue (const std::string& szInput, long& lV);
int GetFloatValue (const std::string& szInput, float& fV);
int GetDoubleValue (const std::string& szInput, double& dV);

// pause-like function
void PAKTC ();

// string manipulations
void TrimLeft (std::string& szTemp);
void TrimRight (std::string& szTemp);
void Trim (std::string& szTemp);
void UpperCase (std::string& szTemp);
void LowerCase (std::string& szTemp);

/* data type limits (see LIMITS.H and FLOATS.H)
 
FLT_MIN             minimum positive value
DBL_MIN             

FLT_MIN_10_EXP      Minimum negative integer such that 10 raised
DBL_MIN_10_EXP      to that number is a representable
                    floating-point number.

FLT_MAX_10_EXP      Maximum integer such that 10 raised to
DBL_MAX_10_EXP      that number is a representable
                    floating-point number.

FLT_MAX             Maximum representable floating-point number.
DBL_MAX

INT_MAX             Maximum (signed) int value 
INT_MIN             Minimum (signed) int value 
LONG_MAX            Maximum (signed) long value 
LONG_MIN            Minimum (signed) long value 

*/

#endif