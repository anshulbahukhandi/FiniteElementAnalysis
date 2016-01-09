/*****************************************************
                A  FEA  Template  Program
            Copyright(c) 2004-11, S. D. Rajan
                  All rights reserved

======================================================
Programming details are in the following document.
======================================================
                       FEAT
         A Finite Element Analysis Template 
      System for Linear Time-Independent Problems


                    S. D. Rajan
          Department of Civil Engineering
              Arizona State University
                   Tempe, AZ 85287
******************************************************/
#include <iostream>
#include "FEAT.h"

int main (int argc, char* argv[])
{
    // the one and only FE object
    CFEAT FEModel; 

    // ------------------------------------- Step 1
    // show program banner
    FEModel.Banner (std::cout);

    // Prepare for I/O
    FEModel.PrepareIO (argc, argv);
    
    // read the problem size
    FEModel.ReadProblemSize ();

    // set problem size
    FEModel.SetSize ();

    // read nodal and element data
    FEModel.ReadModel ();

    // echo the input data in the output file
    FEModel.EchoInput ();

    // ------------------------------------- Steps 2-5
    // construct system equations
    FEModel.ConstructK ();
    FEModel.ConstructF ();

    // ------------------------------------- Step 6
    // impose boundary conditions
    FEModel.ImposeBC ();

    // ------------------------------------- Step 7
    // solve for the nodal displacements
    FEModel.Solve ();

    // ------------------------------------- Steps 8-10
    // compute element response
    FEModel.ComputeResponse ();

    // ------------------------------------- Step 11
    // prints results to output file
    FEModel.CreateOutput ();

    // ------------------------------------- Step 12
    // Close input and output files
    FEModel.TerminateProgram ();

    // inform user
    std::cout << "\nProgram execution completed.\n";

    return 0;
}