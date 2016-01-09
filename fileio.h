/*********************************************
Utility Library Function
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

(non-class) Functions to handle simple file operations

Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_FILEIO_H__
#define __RAJAN_FILEIO_H__

#include <iostream>
#include <fstream>	
#include <string>
#include <vector>

void OpenInputFileByName (const std::string& szPrompt, 
                          std::ifstream& IFile, 
                          const std::ios::openmode&);
void OpenOutputFileByName (const std::string& szPrompt,
                           std::ofstream& OFile,
                           const std::ios::openmode&);
void Rewind (std::ifstream& IOFile);

#endif
