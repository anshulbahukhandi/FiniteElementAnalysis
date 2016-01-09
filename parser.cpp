/*********************************************
Planar Frame Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CParser class implementation.

*********************************************/
#include <cctype>
#include "parser.h"

CParser::CParser ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

CParser::~CParser ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

bool CParser::ReadNextLine (std::ifstream& FileInput, int& nLineNum, 
                            std::string& szInputString, const int MAXCHARS,
                            const std::string& szComment, bool bLowerCase)
// ---------------------------------------------------------------------------
// Function: reads the next line skipping over the comment lines
//           and converts all alphabets to lower case if requested
// Input:    file istream, line #, string to hold the input line,
//           max. # of characters expected in each input line,
//           comment character(s) at the beginning of a comment line,
//           lowercase conversion option
// Output:   updated values of line # and the string
//           return value is true if successful
//                           false if an error state is encountered
// Restriction: Cannot read a line over 256 characters
// ---------------------------------------------------------------------------
{
    int i;
    const int MAXCH = 256;
    static char szInp[MAXCH+1];

    // enough capacity to read and store?
    if (MAXCHARS > MAXCH)
        return false;

    // comment character(s)
    int nCLen = static_cast<int>(szComment.length());

    // read the line (skip over comment lines)
    for(;;)
    {
        ++nLineNum;
        FileInput.getline (szInp, MAXCHARS);

        // end-of-file?
        if (FileInput.eof())
            return false;
        if (FileInput.fail())
            FileInput.clear (FileInput.rdstate() & ~std::ios::failbit);
        // unrecoverable error?
        if (FileInput.bad())
            return false;

        // successful read
        szInputString = szInp;
        if (szInputString.substr(0,nCLen) != szComment)
        {
            // convert to lower case?
            if (bLowerCase)
            {
                for (i=0; i < static_cast<int>(szInputString.length()); i++)
                    szInputString[i] = tolower(szInputString[i]);
            }
            break;
        }
    }

    return true;
}

bool CParser::GetTokens (std::ifstream& IFile,
                         std::vector<std::string>& szVTokens,
                         int& nTokens, const std::string& szDelimiters)
// ---------------------------------------------------------------------------
// Function: reads a line of input from a file and parses the input line
//           into distinct tokens based on the specified delimiter(s)
// Input:    ifstream object, delimiters
// Output:   vector containing the tokens, number of tokens
// ---------------------------------------------------------------------------
{
    const int MAXCHARS = 256;
    static char szEntireLine[MAXCHARS+1];

    IFile.getline (szEntireLine, MAXCHARS);
    if (IFile.eof() || IFile.fail())
        return false;

    // convert input line into std::string object
    std::string szLine (szEntireLine);

    // prepare to store the tokens
    szVTokens.clear();
    std::string::size_type start_loc, end_loc;

    // get location of the first character that is not a delimiter
    start_loc = szLine.find_first_not_of (szDelimiters);

    // loop while the beginning index is not the end of string
    while (start_loc != std::string::npos)
    {
        // get location of the next delimiter
        end_loc = szLine.find_first_of (szDelimiters, start_loc);

        // if this location is the end of string then adjust the value
        // as string length
        if (end_loc == std::string::npos) end_loc = szLine.length();

        // save the string between start_loc and end_loc
        szVTokens.push_back (szLine.substr(start_loc,end_loc-start_loc));

        // get location of the next character that is not a delimiter
        start_loc = szLine.find_first_not_of (szDelimiters, end_loc);
    }
    nTokens = static_cast<int>(szVTokens.size());

    return true;
}
