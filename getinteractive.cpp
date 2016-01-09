/*********************************************
Improved GetInteractive Function
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Object-Oriented Numerical Analysis
*********************************************/
#include "getinteractive.h"
#include <algorithm>
#include <iostream>
#include <cctype>

const int MAXCHARS = 256;	// max. # of characters in input

//=======================================================
//=======================================================
//=================== helper functions ==================
//=======================================================
//=======================================================
// trim leading zeros
void TrimLeadingZeros (std::string& szTemp)
{
    int nPos = 0;		// current parsing position

    // find first nonzero character
    while (szTemp[nPos] == '0' && nPos < MAXCHARS)
        nPos++;

    // length of string
    int nLength = static_cast<int>(szTemp.length());

    // zero-filled string?
    if (nPos == (MAXCHARS-1) || (nLength-nPos) == 0)
    {
        szTemp = "";
        return;
    }

    // extract string without leading zeros
    szTemp = szTemp.substr(nPos, nLength-nPos);
}

// trim leading blanks
void TrimLeft (std::string& szTemp)
{
    int nPos = 0;		// current parsing position

    // find first nonblank character
    while (szTemp[nPos] == ' ' && nPos < MAXCHARS)
        nPos++;

    // length of string
    int nLength = static_cast<int>(szTemp.length());

    // empty string?
    if (nPos == (MAXCHARS-1) || (nLength-nPos) == 0)
    {
        szTemp = "";
        return;
    }

    // extract string without leading blanks
    szTemp = szTemp.substr(nPos, nLength-nPos);

}

// trim trailing blanks
void TrimRight (std::string& szTemp)
{
    int nLast = 0;		// current parsing position

    // find last nonblank character
    int nLength = static_cast<int>(szTemp.length());
    for (int i=nLength-1; i >= 0; i--)
    {
        nLast = i;
        if (szTemp[nLast] != ' ') break;
    }

    // full string?
    if (nLast == (MAXCHARS-1)) {
        return;
    }

    // extract string without trailing blanks
    szTemp = szTemp.substr(0, nLast+1);
}

// trim leading and trailing blanks
void Trim (std::string& szTemp)
{
    TrimLeft (szTemp);
    TrimRight (szTemp);
}

void UpperCase (std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), std::toupper);
}

void LowerCase (std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), std::tolower);
}

int NumChars (long lV, std::string& szT)
{
    int nC = 0;
    long lVT = lV;

    while (lVT > 0)
    {
        szT = szT + ' ';
        lVT /= 10;
        nC++;
    }

    int nX = nC-1;
    while (lV > 0)
    {
        int n = lV - (lV/10)*10;
        lV /= 10;
        szT[nX] = n + '0';
        nX--;
    }

    return nC;
}

// convert alpha string to long integer
int AtoL (const std::string& szTemp, long& lV,
          int nSign)
{
    int nC = 0;
    lV = 0;
    int nLen = static_cast<int>(szTemp.length());

    // capacity of a long number
    //long lVMAX = (nSign >= 1 ? LONG_MAX : labs(LONG_MIN));
    long lVMAX = (nSign >= 1 ? LONG_MAX : LONG_MAX);
    std::string szT;
    int nChars = NumChars(lVMAX, szT);
    if (nLen > nChars)
        return (2);
    else if (nLen == nChars) {
        if (szT < szTemp) return (2);
    }
        
    // Example: value of xyz is
    //          z*10**0 + y*10**1 + z*10**2
    for (int i=nLen-1; i >= 0; i--)
    {
        lV += static_cast<long>((szTemp[i] - '0')*pow(10.0, nC));
        nC++;
    }

    return 0;
}

// convert alpha string to double value < 1
double AtoDFraction (const std::string& szTemp)
{
    double dVF = 0.0;
    int nC = 1;
    // Example: value of .xyz is
    //          x/10**1 + y/10**2 + z/10**3
    for (int i=0; i < static_cast<int>(szTemp.length()); i++)
    {
        dVF += static_cast<double>(szTemp[i] - '0')*pow(10.0, -nC);
        nC++;
    }

    return dVF;
}

// display error message
template <class T>
void ErrorHandler (int nState, T LL=0, T LU=0)
{
    // check the state
    if (nState == 1)
        std::cout << "Error: Invalid input.\n";
    else if (nState == 2)
        std::cout << "Error: Value should be between " << LL
                  << " and " << LU << ".\n";
    else if (nState == 3)
        std::cout << "Error: Check specified range values.\n";
    else if (nState == 4)
        std::cout << "Error: Number of characters restricted to " 
                  << LL << ".\n";
}

// display prompt and read the input
void DisplayandRead (const std::string& szPrompt, 
                     std::string& szInput)
{
    char szRawInput[MAXCHARS + 1];	// actual user input stored here
    int nState;	                    // input state
    std::string szTemp;	            // temporary string

    do {
        // display the prompt
        std::cout << szPrompt;

        // flush any characters in buffer
        std::cin.sync();
        
        // grab all the input until end-of-line
        std::cin.getline(szRawInput, MAXCHARS, '\n');
        szInput = szRawInput;
        nState = std::cin.fail(); // state is non-zero if cin fails
        if (nState)
        {
            std::cin.clear ();
            std::cin >> szTemp;	// since cin failed store the invalid
                                // input in szTemp
        }
    } while (nState);
}

//=======================================================
//=======================================================
//===================== int version =====================
//=======================================================
//=======================================================
void GetInteractive (const std::string& szPrompt,
                     int nUV[], int nSize)
{
    std::string szInput;
    std::string szCurrentString;
    bool bError;
    int nState;

    do 
    {
        DisplayandRead (szPrompt, szInput);
        Trim (szInput);
    
        int i;
        int nPos=0;
        long lV;
        bError = false;
        for (i=0; i < nSize; i++)
        {
            int nLoc = static_cast<int>(szInput.find (" ", nPos));
            szCurrentString = szInput.substr (nPos, nLoc-nPos+1);
            nPos = nLoc+1;
            if ((nState = GetLongValue (szCurrentString, lV)) != 0)
            {
                ErrorHandler (nState, 0, 0);
                bError = true;
                break;
            }
            nUV[i] = static_cast<int>(lV);
        }
    } while (bError);
}

void GetInteractive (const std::string& szPrompt, int& nUV)
{
    GetInteractive (szPrompt, nUV, INT_MIN, INT_MAX);
}

void GetInteractive (const std::string& szPrompt, int& nUV,
                     int nL, int nU)
{
    // range correct?
    if (nL < INT_MIN || nU > INT_MAX)
        ErrorHandler (3, 0, 0);
    
    // initialize
    long lL = nL;
    long lU = nU;
    long lUV;

    // call the long version
    GetInteractive (szPrompt, lUV, lL, lU);

    // convert to int and return value to user
    nUV = static_cast<int>(lUV);
}

//=======================================================
//=======================================================
//==================== int version =====================
//=======================================================
//=======================================================
int GetIntValue (const std::string& szUserInput,
                 int& nV)
{
    long lV;
    int nState;

    if ((nState = GetLongValue (szUserInput, lV)) != 0)
        return 1;
    nV = static_cast<int> (lV);
    return 0;
}

//=======================================================
//=======================================================
//==================== long version =====================
//=======================================================
//=======================================================
int GetLongValue (const std::string& szUserInput,
                  long& lV)
{
    // initialize
    int nSign = 1;	    // assume positive input
    int nPos = 0;	    // current parsing position
    std::string szTemp;	// temporary string
    int nState = 0;     // error state

    // trim leading and trailing blanks
    std::string szInput = szUserInput;
    Trim (szInput);

    // first character + or -
    if (szInput[nPos] == '-') {
        nSign = -1;		// negative integer
        nPos++;
    }
    else if (szInput[nPos] == '+') {
        nPos++;			// ignore + sign
    }

    // extract string without + or -
    int nLast = static_cast<int>(szInput.length());
    szTemp = szInput.substr(nPos, nLast-nPos);
    // trim leading zeros
    if (nLast > 1) TrimLeadingZeros (szTemp);

    // rest of the string is legal?
    int nLoc = static_cast<int>(szTemp.find_first_not_of ("0123456789"));
    if (nLoc >= 0)		// invalid input
        nState = 1;		// set status for invalid input
    // invalid input
    else if (szTemp.length() == 0)
        nState = 1;
    else
    {	
        // legal input so far.
        nState = AtoL (szTemp, lV, nSign);
        // account for proper sign
        lV *= nSign;
    }

    return nState;
}

void GetInteractive (const std::string& szPrompt, long& lUV,
                     long lL, long lU)
{
    // store user input as character string
    std::string szInput;	// parsed user input
    long lV;  		        // store the parsed value here
    int nState;	            // input state

    // valid range?
    if (lL < LONG_MIN || lU > LONG_MAX)
        ErrorHandler (3, 0, 0);
    if (lL >= lU) {
        ErrorHandler (3,0,0);
        return;
    }

    // loop until input is valid
    do {
        nState = 0; // no error

        // display the prompt and read
        DisplayandRead (szPrompt, szInput);

        // get long value
        nState = GetLongValue (szInput, lV); 

        // range check necessary?
        if (nState == 0 && lU > lL) {
            if (lV < lL || lV > lU)
                nState = 2;			// invalid value
        }

        // check the state
        ErrorHandler (nState, lL, lU);

    } while (nState);	// iterate until the input is valid.
                        // state is zero for a valid input

    // return the value to the user
    lUV = lV;
}

void GetInteractive (const std::string& szPrompt, 
                     long& lUV)
{
    GetInteractive (szPrompt, lUV, LONG_MIN, LONG_MAX);
}

//=======================================================
//=======================================================
//=================== float version =====================
//=======================================================
//=======================================================
int GetFloatValue (const std::string& szUserInput,
                   float& fV)
{
    double dV;
    int nState;

    if ((nState = GetDoubleValue (szUserInput, dV)) != 0)
        return 1;
    fV = static_cast<float> (dV);

    return 0;
}

//=======================================================
//=======================================================
//===================== double version ==================
//=======================================================
//=======================================================
void GetInteractive (const std::string& szPrompt, 
                     double& dUV)
{
    GetInteractive (szPrompt, dUV, -DBL_MAX, DBL_MAX);
}

int GetDoubleValue (const std::string& szUserInput, double& dV)
{
    // initialize
    double dSignM = 1.0;	// assume positive mantissa
    double dSignE = 1.0;	// assume positive exponent
    int nPos = 0;		    // current parsing position
    std::string szTemp;     // temporary string
    int nState = 0;         // error state
    // store exponent, mantissa and fractional components
    std::string szExponent= "", szMantissa = "",
                szFraction = "";
    long lVE, lVM; 
    double dVF;

    // trim leading and trailing blanks
    std::string szInput = szUserInput;
    Trim (szInput);

    // first character + or -
    if (szInput[nPos] == '-') {
        dSignM = -1.0;		// negative 
        nPos++;
    }
    else if (szInput[nPos] == '+') {
        nPos++;			// ignore + sign
    }

    // extract string without + or - 
    szTemp = szInput.substr(nPos, szInput.length()-nPos);
    int nLast = static_cast<int>(szTemp.length());
    // trim leading zeros (unless the input is 0)
    if (nLast > 1) TrimLeadingZeros (szTemp);

    // is there an exponent?
    int nLocExp = static_cast<int>(szTemp.find_first_of ("eE"));
    if (nLocExp >= 0) { // yes
        szExponent = szTemp.substr(nLocExp+1, nLast-nLocExp);
        szMantissa = szTemp.substr(0, nLocExp);
        nLast = static_cast<int>(szExponent.length());
        // first character + or - in the exponent
        if (szExponent[0] == '-') {
            dSignE = -1.0;		// negative 
            szExponent = szExponent.substr(1, nLast-1);
        }
        else if (szExponent[0] == '+') {
            szExponent = szExponent.substr(1, nLast-1);
        }
    }
    else // no
        szMantissa = szTemp;

    // is there a decimal?
    int nLocDec = static_cast<int>(szMantissa.find ("."));
    if (nLocDec >= 0) {
        nLast = static_cast<int>(szMantissa.length());
        szFraction = szMantissa.substr(nLocDec+1, nLast-1);
        szMantissa = szMantissa.substr(0, nLocDec);
    }

    // valid string?
    int nLocE = static_cast<int>(szExponent.find_first_not_of ("0123456789"));
    int nLocM = static_cast<int>(szMantissa.find_first_not_of ("0123456789"));
    int nLocF = static_cast<int>(szFraction.find_first_not_of ("0123456789"));
    // invalid input
    if (nLocE >= 0 || nLocM >= 0 || nLocF >= 0)
        nState = 1;		
    // invalid input
    else if (szMantissa.length() == 0 && szFraction.length() == 0)
    {
        dV = 0.0;
        return nState;
    }
    // valid input
    else
    {	
        int nCharsE = static_cast<int>(szExponent.length()); 
        int nCharsM = static_cast<int>(szMantissa.length()); 
        int nCharsF = static_cast<int>(szFraction.length()); 
        // trim leading zeros
        if (nCharsE > 0) {
            TrimLeadingZeros (szExponent);
            nCharsE = static_cast<int>(szExponent.length()); 
        }
        if (nCharsM > 0) {
            TrimLeadingZeros (szMantissa);
            nCharsM = static_cast<int>(szMantissa.length()); 
        }
        // value of exponent
        if ((nState = AtoL (szExponent, lVE, static_cast<int>(dSignE))) != 0)
            return (nState);
        // capacity of a long number
        std::string szT;
        int nChars = NumChars(DBL_MAX_10_EXP, szT);
        if (nCharsE > nChars)
            return (2);
        else if (nCharsE == nChars) {
            if (szT < szExponent) return (2);
        }
            
         // value of mantissa
         if ((nState = AtoL (szMantissa, lVM, static_cast<int>(dSignM))) != 0)
             return (nState);
         // value of fractional component
         dVF = AtoDFraction (szFraction);

         // account for proper sign and final value
         dV = dSignM*(double(lVM) + dVF)*pow(10.0,dSignE*lVE);
    }

    return nState;
}

void GetInteractive (const std::string& szPrompt, double& dUV,
                     double dL, double dU)
{
    // valid input is of the form.
    // [+-][mantissa].[fraction][+-][exponent]
    // [..] is optional

    // store user input as character string
    std::string szInput;  // parsed user input
    double dV;            // store the parsed value here
    int nState = 0;		  // input state

    // valid range?
    if (dL < -DBL_MAX || dU > DBL_MAX)
    {
        ErrorHandler (3,0,0);
        return;
    }
    if (dL >= dU)
    {
        ErrorHandler (3,0,0);
        return;
    }

    // loop until input is valid
    do {
        nState = 0; // no error

        // display the prompt and read
        DisplayandRead (szPrompt, szInput);

        // get double value
        nState = GetDoubleValue (szInput, dV); 

        // range check necessary?
        if (nState == 0 && dU > dL) {
            if (dV < dL || dV > dU)
                    nState = 2;		// invalid value
        }

        // check the state
        ErrorHandler (nState, dL, dU);

    } while (nState);	// iterate until the input is valid.
                        // state is zero for a valid input

    // return the value to the user
    dUV = dV;
}

//=======================================================
//=======================================================
//==================== float version ====================
//=======================================================
//=======================================================
void GetInteractive (const std::string& szPrompt,
                     float& fUV)
{
    GetInteractive (szPrompt, fUV, -FLT_MAX, FLT_MAX);
}

void GetInteractive (const std::string& szPrompt,
                     float& fUV, float fL, float fU)
{
    // call the double version
    double dL = fL;
    double dU = fU;
    double dUV;

    GetInteractive (szPrompt, dUV, dL, dU);

    // convert to float and return value to user
    fUV = static_cast<float>(dUV);
}

//=======================================================
//=======================================================
//==================== string version ===================
//=======================================================
//=======================================================
void GetInteractive (const std::string& szPrompt, float fUV[],
                     int nSize)
{
    std::string szInput;
    std::string szCurrentString;
    bool bError;
    int nState;

    do 
    {
        DisplayandRead (szPrompt, szInput);
        Trim (szInput);
    
        int i;
        int nPos=0;
        double dV;
        bError = false;
        for (i=0; i < nSize; i++)
        {
            int nLoc = static_cast<int>(szInput.find (" ", nPos));
            szCurrentString = szInput.substr (nPos, nLoc-nPos+1);
            nPos = nLoc+1;
            if ((nState = GetDoubleValue (szCurrentString, dV)) != 0)
            {
                ErrorHandler (nState, 0, 0);
                bError = true;
                break;
            }
            fUV[i] = static_cast<float>(dV);
        }
    } while (bError);
}

void GetInteractive (const std::string& szPrompt, std::string& Value,
                     int nMaxChars)
{
    int nState;
    std::string szTemp;
    char szRawInput[MAXCHARS + 1];	// actual user input stored here

    do {
        std::cout << szPrompt;

        // flush any characters in buffer
        std::cin.sync();
        
        // grab all the input until end-of-line
        std::cin.getline (szRawInput, MAXCHARS, '\n');
        nState = std::cin.fail(); // state is non-zero if cin fails
        if (nState)
        {
            std::cin.clear ();
            std::cin >> szTemp;	// since cin failed store the invalid
                               // input in szTemp
        }
        else
        {
            Value = szRawInput;
            if (static_cast<int>(Value.length()) > nMaxChars) 
                nState = 4;
        }

        // check the state
        ErrorHandler (nState, nMaxChars, 0);

    } while (nState);	// iterate until the input is valid.
                        // state is zero for a valid input
}

void FlushInput ()
{
    std::string szTemp;

    std::cin.clear (0);
    std::cin >> szTemp;	// since cin failed store the invalid
}

// pause-like function
void PAKTC ()
{
    std::cout << "Press any key to continue.";
    char ch;
    std::cin.get (ch);

    // flush the input buffer just in case
    while (!std::cin.eof() && ch != '\n')
        std::cin.get (ch);
}

