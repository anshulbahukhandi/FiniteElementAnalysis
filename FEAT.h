/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-06, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_CFEAT_H__
#define __RAJAN_CFEAT_H__

#include <fstream>              // needed for file input and output
#include <string>               // std::string class

#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"  // vector container
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\matrixtemplate.h"  // matrix container
#include "constants.h"                  // program limitations etc.
#include "node.h"                       // node class
#include "element.h"                    // element class
#include "nodalresponse.h"              // nodal response class
#include "elementresponse.h"            // element response class
#include "elementloads.h"

class CFEAT
{
    public:
        CFEAT ();
        ~CFEAT ();

        enum ErrorStates {NUMNODES, NUMELEMENTS, DEBUGCODE,
                          NODENUMBER, ELEMENTNUMBER, EPROPERTY,
                          YOUNGSMODULUS, UNSTABLE, INVALIDINPUT,
                          INVALIDLASTLINE, NODALFIXITY,
                          INPUTFILE, OUTPUTFILE, COMMANDLINE,ELEMENTLOAD,INVERSEERROR};

        // helper functions
        void Banner (std::ostream& OF) const;
        void PrepareIO (int argc, char* argv[]);
        void ReadProblemSize ();
        void ReadModel ();
        void EchoInput ();
        void ConstructK ();
        void ConstructF ();
        void ImposeBC ();
        void Solve ();
        void ComputeResponse ();
        void CreateOutput ();
        void TerminateProgram ();

        // modifier functions
        void SetSize ();

    private:
        std::string m_szDescription;    // problem description
        int m_nNodes;                   // number of nodes
        int m_nElements;                // number of elements
        int m_nDOF;                     // total degrees-of-freedom
        int m_nDebugLevel;              // debugging level
        int m_nLineNumber;              // current line number in input 
		int m_nElementLoads;
		int m_ElementLoad;
		int m_nNodalFixity;
		int m_nNodalloads;

        // these store the FE model data
        CVector<CNode>            m_NodalData;           // nodal data
        CVector<CElement>         m_ElementData;         // element data
        CVector<CNodalResponse>   m_NodalResponseData;   // nodal response
        CVector<CElementResponse> m_ElementResponseData; // element response
		//CVector<CElementLoads>	  m_nElementLoads;
		CVector<CElementLoads>	  m_nElementLoad;
        std::ifstream m_FileInput;  // File Input
        std::ofstream m_FileOutput; // File Output

        CMatrix<double> m_SSM,m_SM;      // structural stiffness matrix
        CMatrix<double> m_SND,m_SD;      // structural nodal displacements
        CMatrix<double> m_SNF,m_SF;      // structural nodal forces

        void ErrorHandler (ErrorStates); // handles detected errors
        void SuppressDOF  (int);         // imposes EBC

        // element-related
        int k1DC0LElement (int nE, CVector<int>& nVNList,
                           CMatrix<double>& dMk,CMatrix<double>& dTF);
        void s1DC0LElement (int nE);
};

#endif