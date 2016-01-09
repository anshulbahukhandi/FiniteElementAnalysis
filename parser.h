/*********************************************
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CParser class definition.

*********************************************/
#ifndef __RAJAN_PARSER_H__
#define __RAJAN_PARSER_H__

#include <fstream>
#include <string>
#include <vector>

class CParser
{
    public:
        CParser ();
        ~CParser ();

        bool ReadNextLine (std::ifstream& FileInput, int& nL,
                           std::string& szInputString,
                           const int MAXCHARS, 
                           const std::string& szComment,
                           bool bLowerCase = true);
        bool GetTokens (std::ifstream& IFile,
                        std::vector<std::string>& szVTokens,
                        int& nTokens, const std::string& szDelimiters);

    private:

};

#endif