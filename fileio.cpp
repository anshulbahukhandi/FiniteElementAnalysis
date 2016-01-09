/*********************************************
Utility Function
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

(non-class) Functions to handle simple file operations

Object-Oriented Numerical Analysis
*********************************************/
#include "fileio.h"
#include "getinteractive.h"

void OpenInputFileByName (const std::string& szPrompt,
                          std::ifstream& IFile,
                          const std::ios::openmode& oMode)
// ---------------------------------------------------------------------------
// Function: opens a file for reading
// Input:    prompt, ifstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
	bool bDone = false;
    std::string szFileName;

	do
    {
        GetInteractive (szPrompt, szFileName, 256);            //!!!!what is this szfilename and szprompt for??

		// open the file
		IFile.open (szFileName.c_str(), oMode); 
		if (!IFile)
        {
            std::cout << "Unable to open file." << std::endl;
			IFile.clear ();
		}
		else
			bDone = true; // file opened successfully
	} while (!bDone);
}

void OpenOutputFileByName (const std::string& szPrompt,
                           std::ofstream& OFile,
                           const std::ios::openmode& oMode)
// ---------------------------------------------------------------------------
// Function: opens a file for writing
// Input:    prompt, ofstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
	bool bDone = false;
    std::string szFileName;

	do
    {
        GetInteractive (szPrompt, szFileName, 256);

		// open the file
		OFile.open (szFileName.c_str(), oMode); 
		if (!OFile)
        {
			std::cout << "Unable to open file." << std::endl;
			OFile.clear ();
		}
		else
			bDone = true; // file opened successfully
	} while (!bDone);
}

void Rewind (std::ifstream& IOFile)                                                  //what is this for?????
// ---------------------------------------------------------------------------
// Function: rewinds an input file
// Input:    ifstream object
// Output:   file positioned at the beginning 
// ---------------------------------------------------------------------------
{
    // clear all error bits
    IOFile.clear (std::ios_base::goodbit);

    // now position the file at byte 0 (beginning of the file)
    IOFile.seekg (0L, std::ios::beg);
}

