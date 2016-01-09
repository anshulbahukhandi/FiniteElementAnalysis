/***********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved

------------------------------------------------
This version is for a 2D solid mechanics problem
------------------------------------------------

***********************************************/
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "FEAT.h"
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\fileio.h"
#include "F:\ASU Courses\cee532\Library\Library\MatToolBox.h"
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\parser.h"
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\getinteractive.h"

std::string szInputString;
std::string szDelimiters(",");
std::string szComment = "**";


CFEAT::CFEAT ()
// ==================================================================
// Function: Constructor
//    Input: None
//   Output: None
// ==================================================================
{
    m_nNodes      = 0;
    m_nElements   = 0;
    m_nDOF        = 0;
    m_nDebugLevel = 0;
    m_nLineNumber = 0;
	m_ElementLoad = 0;
	m_nNodalFixity = 0;
	m_nNodalloads = 0;
}

CFEAT::~CFEAT ()
// ==================================================================
// Function: Destructor
//    Input: None
//   Output: None
// ==================================================================
{
}

void CFEAT::Banner (std::ostream& OF) const
// ==================================================================
// Function: Prints the banner on output stream
//    Input: Output stream
//   Output: None
// ==================================================================
{
    OF << "\t\t--------------------------------------------\n";
    OF << "\t\t           A  FEA  Template  Program        \n";
    OF << "\t\t         Finite Elements for Engineers      \n";
    OF << "\t\t           (c) 2000-11, S. D. Rajan         \n";
    OF << "\t\t          Enhanced By: Anshul Bahukhandi      \n";
    OF << "\t\t--------------------------------------------\n";
}

void CFEAT::PrepareIO (int argc, char* argv[])
// ==================================================================
// Function: Opens the input and output files
//    Input: Number of command line arguments, list of arguments
//   Output: None
// ==================================================================
{
    // no command line arguments. hence get file names interactively
    if (argc <= 1)
    {
        // open the input file
        OpenInputFileByName ("Complete input file name: ", 
                             m_FileInput, std::ios::in);

        // open the output file
        OpenOutputFileByName ("Complete output file name: ", 
                              m_FileOutput, std::ios::out);
    }
    // user has specified command line arguments
    else if (argc == 3)
    {
        m_FileInput.open (argv[1], std::ios::in); 
		if (!m_FileInput)
			ErrorHandler (INPUTFILE);
        m_FileOutput.open (argv[2], std::ios::out); 
		if (!m_FileOutput)
			ErrorHandler (OUTPUTFILE);
        std::cout << "\nReading FE model from " << argv[1] << '\n';
        std::cout << "\nCreating results in   " << argv[2] << '\n';
    }
    else
    {
        ErrorHandler (COMMANDLINE);
    }

    // print banner
    Banner (m_FileOutput);
}

void CFEAT::ReadProblemSize ()
// ==================================================================
// Function: Reads the problem size from the input file
//    Input: None
//   Output: None
// ==================================================================
{
    CParser Parse;

    // read the problem description
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (INVALIDINPUT);
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (INVALIDINPUT);

    // nodal coordinates
    if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                             MAXCHARS, szComment))
        ErrorHandler (INVALIDINPUT);
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            ErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,13) == "*nodal fixity")
            break;
        ++m_nNodes;
    }

    // nodal fixity
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            ErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,12) == "*nodal loads")
            break;
		++m_nNodalFixity;
    }

	 for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            ErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,13) == "*element data")
            break;
		++m_nNodalloads;
    }
    // element data
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            ErrorHandler (INVALIDINPUT);
        if (szInputString.substr(0,14) == "*element loads")
            break;
        ++m_nElements;
    }

    // element loads
    for (;;)
    {
        if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
                                 MAXCHARS, szComment))
            ErrorHandler (INVALIDINPUT);
            if (szInputString.substr(0,4) == "*end")
                break;
        ++m_ElementLoad;
    }

    // check data for validity
    if (m_nNodes <= 1) 
        ErrorHandler (NUMNODES);
    if (m_nElements <= 0) 
        ErrorHandler (NUMELEMENTS);
    if (m_nDebugLevel < 0 || m_nDebugLevel > 1) 
        ErrorHandler (DEBUGCODE);
}

void CFEAT::SetSize ()
// ==================================================================
// Function: Allocates memory for all the major vectors and matrices
//           used in the program
//    Input: None
//   Output: None
// ==================================================================
{
    // allocate space for nodal data
    m_NodalData.SetName ("NodalData");
    m_NodalData.SetSize (m_nNodes);

    // allocate space for nodal response data
    m_NodalResponseData.SetName ("NodalResponseData");
    m_NodalResponseData.SetSize (m_nNodes);

    // allocate space for element data
    m_ElementData.SetName ("ElementData");
    m_ElementData.SetSize (m_nElements);

    // allocate space for element response data
    m_ElementResponseData.SetName ("ElementResponseData");
    m_ElementResponseData.SetSize (m_nElements);
	if(m_ElementLoad>0)
	{
	m_nElementLoad.SetName ("ElementLoads");
    m_nElementLoad.SetSize (m_ElementLoad);
	}
    // allocate and initialize major matrices
    m_nDOF = DOFPN*m_nNodes;     // total number of unknowns
    m_SSM.SetName ("SSM");       // system (structural stiffness) matrix
    m_SSM.SetSize (m_nDOF, m_nDOF);     m_SSM.Set (0.0);
	m_SM.SetSize (m_nDOF, m_nDOF);     m_SM.Set (0.0);
    m_SNF.SetName ("SNF");       // nodal forces
    m_SNF.SetSize (m_nDOF, 1);          m_SNF.Set (0.0);
	 m_SF.SetSize (m_nDOF, 1);          m_SF.Set (0.0);
    m_SND.SetName ("SND");       // nodal unknowns
    m_SND.SetSize (m_nDOF, 1);          m_SND.Set (0.0);
	 m_SD.SetSize (m_nDOF, 1);          m_SD.Set (0.0);
}

void CFEAT::ReadModel ()
// ==================================================================
// Function: Reads the FE model details
//    Input: None
//   Output: None
// ==================================================================
{
	 CParser Parse;
    std::vector<std::string> szVTokens;
    int nTokens;
    int i, j,nTag;
	Rewind (m_FileInput);
	if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
        ErrorHandler (INVALIDINPUT);
    if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
        ErrorHandler (INVALIDINPUT);
    // *nodal coordinates
    if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
        ErrorHandler (INVALIDINPUT);
    CVector<float> fVCoor("fVCoor", PDIM);
    for (i=1; i <= m_nNodes; i++)
    {
        if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
            ErrorHandler (INVALIDINPUT);
        if (GetIntValue (szVTokens[0], nTag) != 0)
            ErrorHandler (INVALIDINPUT);
        if (nTag <= 0 || nTag > m_nNodes)
            ErrorHandler (NODENUMBER);
        for (j=1; j <= PDIM; j++)
        {
            if (GetFloatValue (szVTokens[j], fVCoor(j)) != 0)
                ErrorHandler (INVALIDINPUT);
        }

        // store the data
        m_NodalData(nTag).SetCoords (fVCoor);
    }

    // *nodal fixity
    if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
        ErrorHandler (INVALIDINPUT);
    CVector<CNode::NBConditions> VFC(DOFPN);
    CVector<float>VDisp(DOFPN);
	std::string szXFC,szYFC;
	float fXDisp,fYDisp;
    for (i=1; i <= m_nNodalFixity; i++)
    {
        if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
            ErrorHandler (INVALIDINPUT);
        if (GetIntValue (szVTokens[0], nTag) != 0)
            ErrorHandler (INVALIDINPUT);
        if (nTag <= 0 || nTag > m_nNodes)
            ErrorHandler (NODENUMBER);
        // fixity conditions
		if(szVTokens[1]==" known")
			szXFC = "known";
		else if(szVTokens[1]==" free")
			szXFC = "free";
		else
            ErrorHandler (NODALFIXITY);	
		if(szVTokens[2]==" known")
			szYFC = "known";
		else if(szVTokens[2]==" free")
			szYFC = "free";
		else
            ErrorHandler (NODALFIXITY);
			if (GetFloatValue (szVTokens[3], fXDisp) != 0)
                ErrorHandler (INVALIDINPUT);
			if (GetFloatValue (szVTokens[4], fYDisp) != 0)
                ErrorHandler (INVALIDINPUT);
        VFC(1) = (szXFC == "known" ? CNode::KNOWN : CNode::FREE);
		VFC(2) = (szYFC == "known" ? CNode::KNOWN: CNode::FREE);
		VDisp(1)=fXDisp;VDisp(2)=fYDisp;
        m_NodalData(nTag).SetFixity (VFC);
		m_NodalData(nTag).SetDisps (VDisp);	        
    }
	//nodal loads
	if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
        ErrorHandler (INVALIDINPUT);
	CVector<float> fVForces(DOFPN);
	float fVtemp;
		 for (i=1; i <= m_nNodalloads; i++)
		 {
			 if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
				ErrorHandler (INVALIDINPUT);
			if (GetIntValue (szVTokens[0], nTag) != 0)
			    ErrorHandler (INVALIDINPUT);
			if (nTag <= 0 || nTag > m_nNodes)
            ErrorHandler (NODENUMBER);
			   for (j=1; j <= DOFPN; j++)
				{
				if (GetFloatValue (szVTokens[j], fVForces(j)) != 0)
					ErrorHandler (INVALIDINPUT);
				}
			   if (GetFloatValue (szVTokens[DOFPN+1], fVtemp) != 0)
					ErrorHandler (INVALIDINPUT);
        int nSDOF = (nTag-1)*DOFPN;
        for (int l=1; l <=DOFPN; l++)
            m_SNF(nSDOF+l,1) = fVForces(l);
		m_NodalData(nTag).SetLoads(fVForces);
		m_NodalData(nTag).SetTemp(fVtemp);
		 }
    // *element data
    if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
        ErrorHandler (INVALIDINPUT);
    CVector<float> fVMatP(NUMEP);
	CElement::type VC;
	std::string szStr;
	CVector <int> nVNList;
	nVNList.SetSize (NUMENODES);
    for (i=1; i <= m_nElements; i++)
    {
        // list of nodes
        if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
            ErrorHandler (INVALIDINPUT);
		if (GetIntValue (szVTokens[0],nTag) != 0)
            ErrorHandler (INVALIDINPUT);
		if (nTag <= 0 || nTag > m_nElements)
            ErrorHandler (ELEMENTNUMBER);
		
		if (szVTokens[1]==" plane strain")
			szStr="plane strain";
		else if(szVTokens[1]==" plane stress")
			szStr="plane stress";
		else
			ErrorHandler(INVALIDINPUT);
        for (j=1; j <= NUMENODES; j++)
        {
            if (GetIntValue (szVTokens[j+1], nVNList(j)) != 0)
                ErrorHandler (INVALIDINPUT);
            if (nVNList(j) <= 0 || nVNList(j) > m_nNodes)
                ErrorHandler (NODENUMBER);
        }
		for (int l=1; l <=NUMEP; l++)
        {
            if (GetFloatValue (szVTokens[NUMENODES+l+1],fVMatP(l)) != 0)
                ErrorHandler (INVALIDINPUT);
        }
		if (fVMatP(1) <= 0.0f)
            ErrorHandler (YOUNGSMODULUS);
        // store the data
        VC = (szStr == "plane stress" ? CElement::PLANESTRESS : CElement::PLANESTRAIN);
		m_ElementData(nTag).Setstressstrain (VC);
        m_ElementData(nTag).SetMaterialProperty (fVMatP);
        m_ElementData(nTag).SetNodes (nVNList);
    }

    // *element properties
    if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
            ErrorHandler (INVALIDINPUT);
	float fphase,fVELoad;
    for (i=1; i <= m_ElementLoad; i++)
    {
        // thickness and material properties
        if (!Parse.GetTokens (m_FileInput, szVTokens, nTokens, szDelimiters))
            ErrorHandler (INVALIDINPUT);
        if (GetIntValue (szVTokens[0], nTag) != 0)
            ErrorHandler (INVALIDINPUT);
		if (nTag <= 0 || nTag > m_nElements)
            ErrorHandler (NODENUMBER);
        if (GetFloatValue (szVTokens[1], fphase) != 0)
            ErrorHandler (INVALIDINPUT);
        if (GetFloatValue (szVTokens[2], fVELoad) != 0)
            ErrorHandler (INVALIDINPUT);
		m_nElementLoad(i).SetElementalload (nTag,fphase,fVELoad);
    }
 
}

void CFEAT::ImposeBC ()
// ==================================================================
// Function: Imposes the homogenous EBC for all nodes
//    Input: None
//   Output: None
// ==================================================================
{
    int i, j;
    CVector<CNode::NBConditions> nVNBC(DOFPN);
    
    // loop thro' all nodes
    for (i=1; i <= m_nNodes; i++)
    {
        m_NodalData(i).GetFixity (nVNBC);
        for (j=1; j <= DOFPN; j++)
        {
            // if condition is HOMOGENOUS then modify equation
            if (nVNBC(j) == CNode::KNOWN)
            {
                int nGDOF = DOFPN*(i-1)+j;
                SuppressDOF (nGDOF);
            }
        }
    }
}

void CFEAT::SuppressDOF (int nEqn)
// ==================================================================
// Function: Imposes the homogenous EBC for the specified equation
//    Input: Equation number
//   Output: None
// ==================================================================
{
  
	for (int j=1; j <= m_nDOF; j++)
    {
		m_SNF(j,1)-=(m_SSM(j,nEqn)*m_SND(nEqn,1));
        // zero out the row
        m_SSM(nEqn, j) = 0.0;
        // zero out the column
        m_SSM(j, nEqn) = 0.0;
    }

    m_SSM(nEqn, nEqn) = 1.0;   // set diagonal to 1
    m_SNF(nEqn, 1) = 0.0;      // set RHS to zero
}

void CFEAT::Solve ()
// ==================================================================
// Function: Solves the system equilibrium equations
//    Input: None
//   Output: None
// ==================================================================
{
    int i, j;

    // implement Gauss Elimination Technique
	//if (AxEqb (m_SSM, m_SND, m_SNF, TOL) != 0)
 //       // error because of user input?
 //       ErrorHandler (UNSTABLE);
		 bool nError = AxEqb(m_SSM,m_SND,m_SNF,TOL);
		 if (!nError)
    {
        ErrorHandler (UNSTABLE);
    }
    else
    {
        // no error. store the nodal displacements
        CVector<float> fVResp(DOFPN);
        for (i=1; i <= m_nNodes; i++)
        {
            for (j=1; j <= DOFPN; j++)
            {
                int nGDOF = DOFPN*(i-1)+j;
                fVResp(j) = static_cast<float>(m_SND(nGDOF,1));
            }
            m_NodalResponseData(i).SetValues (fVResp);
        }
    }
		 m_SD=m_SND;
		 Multiply(m_SM,m_SD,m_SF);
}

void CFEAT::ConstructK ()
// ==================================================================
// Function: Generates k and assembles them into K
//    Input: None
//   Output: None
// ==================================================================
{
    int i, j, k;
    const int KSIZE = NUMENODES*DOFPN;     // size of the stiffness matrix
    CVector<int>    nVEDOF(KSIZE);         // dof associated with element
    CMatrix<double> dMk("k",KSIZE,KSIZE);  // to store the element stiffness 
                                           // matrix
	 CMatrix<double> dTF("t",KSIZE,1);
	 dMk.Set(0.0);
	 dTF.Set(0.0);
    CVector<int> nVNL(NUMENODES); // list of element nodes
	
    // loop thro' all elements
    for (i=1; i <= m_nElements; i++)
    {
		m_ElementData(i).GetNodes(nVNL);
        // Generate k for 1D-C0 linear element
        k1DC0LElement (i, nVNL, dMk,dTF);

        // get global degrees-of-freedom associated with element
        int nIndex = 1;
        for (j=1; j <= NUMENODES; j++)
        {
            int n = (nVNL(j)*DOFPN)-1;
            for (k=1; k <= DOFPN; k++)
			{
                nVEDOF(nIndex) = n;
				n++;nIndex++;
				
			}
        }

        // assemble into structural K
        for (j=1; j <= KSIZE; j++)
        {
            int nRow = nVEDOF(j);
            for (k=1; k <= KSIZE; k++)
            {
                int nCol = nVEDOF(k);
                m_SSM(nRow, nCol) += dMk(j,k);
            }
			m_SNF(nRow, 1) += dTF(j,1);
        }
    }
	m_SM=m_SSM;
	
    // debug?
    if (m_nDebugLevel == 1)
    {
        PrintMatrixRowWise (m_SSM, 
            "Structural Stiffness (Before BCs)", m_FileOutput);
    }
}

void CFEAT::ConstructF ()
// ==================================================================
// Function: Constructs the system nodal vector
//    Input: None
//   Output: None
// ==================================================================
{
    int i, E;
    float fv1,fv2,flength;
	CVector <float> fVLoads(2);
	CVector <int>fcord(NUMEP);
	CVector <float> fpt1(PDIM),fpt2(PDIM);

	for (i=1; i <= m_ElementLoad; i++)
    {
		m_nElementLoad(i).GetElementalload (E,fv1,fv2);
		m_ElementData(E).GetNodes(fcord);
		if(fv1==1)
		{
			m_NodalData(fcord(1)).GetCoords(fpt1);
			m_NodalData(fcord(2)).GetCoords(fpt2);
			flength=sqrt(pow((fpt2(1)-fpt1(1)),2) + ((pow (fpt2(2)-fpt1(2),2))));
				m_SNF(2*fcord(1)-1,1)+=((flength*((7*fv2)+(3*fv2)))/20);
				m_SNF(2*fcord(1),1)+=((pow(flength,2)*(3*fv2+2*fv2))/60);
				m_SNF(2*fcord(2)-1,1)+=((flength*((3*fv2)+(7*fv2)))/20);
				m_SNF(2*fcord(2),1)+=-((pow(flength,2)*(2*fv2+3*fv2))/60);
		}
      if(fv1==2)
		{
			m_NodalData(fcord(2)).GetCoords(fpt1);
			m_NodalData(fcord(3)).GetCoords(fpt2);
			flength=sqrt(pow((fpt2(1)-fpt1(1)),2) + ((pow (fpt2(2)-fpt1(2),2))));
				m_SNF(2*fcord(2)-1,1)+=((flength*((7*fv2)+(3*fv2)))/20);
				m_SNF(2*fcord(2),1)+=((pow(flength,2)*(3*fv2+2*fv2))/60);
				m_SNF(2*fcord(3)-1,1)+=((flength*((3*fv2)+(7*fv2)))/20);
				m_SNF(2*fcord(3),1)+=-((pow(flength,2)*(2*fv2+3*fv2))/60);
		}
	  if(fv1==3)
		{
			m_NodalData(fcord(3)).GetCoords(fpt1);
			m_NodalData(fcord(4)).GetCoords(fpt2);
			flength=sqrt(pow((fpt2(1)-fpt1(1)),2) + ((pow (fpt2(2)-fpt1(2),2))));
				m_SNF(2*fcord(3)-1,1)+=((flength*((7*fv2)+(3*fv2)))/20);
				m_SNF(2*fcord(3),1)+=((pow(flength,2)*(3*fv2+2*fv2))/60);
				m_SNF(2*fcord(4)-1,1)+=((flength*((3*fv2)+(7*fv2)))/20);
				m_SNF(2*fcord(4),1)+=-((pow(flength,2)*(2*fv2+3*fv2))/60);
		}
	  if(fv1==4)
		{
			m_NodalData(fcord(4)).GetCoords(fpt1);
			m_NodalData(fcord(1)).GetCoords(fpt2);
			flength=sqrt(pow((fpt2(1)-fpt1(1)),2) + ((pow (fpt2(2)-fpt1(2),2))));
				m_SNF(2*fcord(4)-1,1)+=((flength*((7*fv2)+(3*fv2)))/20);
				m_SNF(2*fcord(4),1)+=((pow(flength,2)*(3*fv2+2*fv2))/60);
				m_SNF(2*fcord(1)-1,1)+=((flength*((3*fv2)+(7*fv2)))/20);
				m_SNF(2*fcord(1),1)+=-((pow(flength,2)*(2*fv2+3*fv2))/60);
		}
	  
    }
	/*std::ostringstream szPrompt;
        szPrompt << "SSF";
		PrintMatrixRowWise (m_SNF, szPrompt.str(), m_FileOutput);
		m_FileOutput << '\n';*/
}

void CFEAT::ComputeResponse ()
// ==================================================================
// Function: Computes the element strain and stress values
//    Input: None
//   Output: None
// ==================================================================
{
    int i;
   
    // loop thro' all elements
    for (i=1; i <= m_nElements; i++)
    {
        // Generate k for 1D-C0 linear element
        s1DC0LElement (i);
    }
}

void CFEAT::EchoInput ()
// ==================================================================
// Function: Creates the contents of the output file
//           echoing the input data
//    Input: None
//   Output: None
// ==================================================================
{
    
    std::string szNDOF[] = {"Known", "Free"};
  int i, j;

	
	CElement::type VC;
	std::string szS;
	m_ElementData(1).Getstressstrain (VC);
	if (VC==0)
		szS="Plane Stress";
	else
		szS="Plane Strain";
	// print the problem description
	m_FileOutput << "                 PROBLEM DESCRIPTION" << '\n';
    m_FileOutput << "                 ===================" << '\n';
	m_FileOutput << "                         JOBNAME : FEM Analysis of Q4 Element"  << '\n';
    m_FileOutput << "              NUMBER OF ELEMENTS = " << m_nElements << '\n';
    m_FileOutput << "                 NUMBER OF NODES = " << m_nNodes << '\n';
	m_FileOutput << "  NUMBER OF ELEMENTAL LOAD CASES = " <<m_ElementLoad << "\n";
	m_FileOutput << "      NUMBER OF NODAL LOAD CASES = " <<m_nNodalloads << "\n";
	m_FileOutput << "					  PROBLEM TYPE= " <<szS << "\n";
	
	CVector<float > fVC1(PDIM);
	CVector<CNode::NBConditions> VFC(DOFPN);
	// print the nodal data
	m_FileOutput << "\n\n";
    m_FileOutput << "NODAL INFORMATION\n";
    m_FileOutput << "NODE"<< "       "<< "X COOR"<< "       "<< "Y COOR"<< "\n";
	m_FileOutput << "--------------------------------------------------------------------------------------------------------\n";
    m_FileOutput << "    "<< "       "<< "  m   "<< "       "<< "  m   "<< "\n";
	m_FileOutput << "--------------------------------------------------------------------------------------------------------\n";
	 for(i=1;i<=m_nNodes;i++)
	 {
		 
		  m_NodalData(i).GetCoords (fVC1);
		  m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << i << "       "<< std::setw(6) << fVC1(1)
					   << "       "<< std::setw(6) << fVC1(2)<<"\n";
	 }
	 m_FileOutput << "\n\n";
    m_FileOutput << "NODAL INFORMATION\n";
    m_FileOutput << "NODE"<< "       "<< "X-FIXITY  "<< "       "<< "Y-FIXITY  " << "\n";
	m_FileOutput << "    "<< "       "<< "          "<< "       "<< "          " << "\n";
	 m_FileOutput << "--------------------------------------------------------------------------------------------------------\n";
	 for(i=1;i<=m_nNodes;i++)
	 {
	 m_NodalData(i).GetFixity (VFC);
	
	  m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << i << "       ";
		 if((VFC(1)==0)&&(VFC(2)==0))
				m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) <<"knonwn"<<  std::setw(10) << "knonwn  \n";
		 else if((VFC(1)==1)&&(VFC(2)==0))
			 m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) << "free"<<  std::setw(10) << "knonwn \n";
		 else if((VFC(1)==0)&&(VFC(2)==1))
			 m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) << "knonwn"<<  std::setw(10) << "free \n";
		 else if((VFC(1)==1)&&(VFC(2)==1))
			  m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) << "free"<<  std::setw(10) << "free \n";
	
	
	 }
	 CVector<float>fVForces(PDIM),VDisp(PDIM);
	 float fTemp;
	 //print the nodal loads and temperature information
	 m_FileOutput << "\n\n";					  
	 m_FileOutput << "NODAL LOADS AND TEMPERATURE CHANGE\n";
	 m_FileOutput << "Node#"<< "       "<<"X-FORCE   "<< "       "<< "Y-FORCE   "<< "       "<< "Temp Change "<< "X-DISP     "<< "Y-DISP     "<< "\n";
	 m_FileOutput << "----------------------------------------------------------------------------------------------------------------------------------\n";
	 m_FileOutput << "    "<< "       "<<"   N      "<< "       "<< "    N     "<< "       "<< "   C      "<< "  m       "<< "  m       "<< "\n";
	 m_FileOutput << "----------------------------------------------------------------------------------------------------------------------------------\n";					  
	 for(i=1;i<=m_nNodes;i++)
	 {
		m_NodalData(i).GetLoads(fVForces);
		m_NodalData(i).GetTemp(fTemp);
		m_NodalData(i).GetDisps (VDisp);	 
		 m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << i << "       ";
		 for(j=1;j<=DOFPN;j++)
			  m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << fVForces(j) << "       ";	
		 m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << fTemp;
		 for(j=1;j<=DOFPN;j++)
			  m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << VDisp(j) << "       ";	
		 m_FileOutput <<"\n";
	 }				  

     //print the elemental loads
	 float fphase,fVELoad;
	 int nE;
	 m_FileOutput << "\n\n";
	 m_FileOutput << "ELEMENTAL LOADS                       \n";
	 m_FileOutput << "ELEM"<< "       "<< "LOAD SIDE   "
							<<"       "<< "LOAD INTENSITY"<< "\n";
	 m_FileOutput << "---------------------------------------------------------------------------------------------\n";		
	 m_FileOutput << "    "<< "       "<< "            "
		                    <<"       "<< "              "<< "\n";
	 m_FileOutput << "    "<<"        "<< "            "
						    <<"       "<< "(N)or(N -m )  "<<"\n";
	m_FileOutput << "---------------------------------------------------------------------------------------------\n";	
	for(i=1;i<=m_nElementLoads;i++)
	{
		std::string Type;
		m_nElementLoad(i).SetElementalload (nE,fphase,fVELoad);
	
		 m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << nE << "       "
					  << std::setw(12) << fphase << "       "<< std::setw(20) << fVELoad<< "\n";
	}	


    //print the Element nodes and properties
	CVector<float> fVMatP(NUMEP);
	
	std::string szStr;
	CVector <int> nVNList(NUMENODES);
	m_FileOutput << "\n\n";
	m_FileOutput << "ELEMENT NODES AND PROPERTIES INFORMATION											 \n";
	m_FileOutput << "--------------------------------------------------------------------------------------------------------\n";
	m_FileOutput << " ELEMENT "<< "       "<< "NODE 1" <<"      " << "NODE 2" <<"      "<< "NODE 3" <<"      "<< "NODE 4" <<"      "<< "    THICKNESS"<< "       "<< "YOUNG'S MODULUS"<< "       "<< "POISSON'S RATIO"
						   << "       "<< "Thermal Coefficient"<<"\n";				
	m_FileOutput << "         "<< "       "<< "       m     "<< "       "<< "     N/m^2     "<< "       "<< "               "
						   << "       "<< "   m/m - C    "<<"\n";
	m_FileOutput << "--------------------------------------------------------------------------------------------------------\n";
	for(i=1;i<=m_nElements;i++)
	{
		
        m_ElementData(i).GetMaterialProperty (fVMatP);
        m_ElementData(i).GetNodes (nVNList);
		m_FileOutput<< setiosflags(std::ios::left) << std::setw(5) << i ;
		for(j=1;j<=NUMENODES;j++)
		m_FileOutput << setiosflags(std::ios::left) <<"       "
					 << std::setw(15) << nVNList(j) <<"       ";
		for(j=1;j<=NUMEP;j++)
		m_FileOutput << setiosflags(std::ios::left) <<"       "
					 << std::setw(15) << fVMatP(j) <<"       ";
		m_FileOutput<<"\n";
	}
	}
   

void CFEAT::CreateOutput ()
// ==================================================================
// Function: Creates the contents of the output file
//           with the results from FE analysis
//    Input: None
//   Output: None
// ==================================================================
{	
	
	int i, j;
	CVector<CNode::NBConditions> VFC(DOFPN);
	CVector<float> fVStress(NUMSCOMP),fVStrain(NUMSCOMP);
	float fS11,fS22,fVonMis,fTersac;
	CVector<float> fVResp(PDIM);
	// print the vital statistics
	m_FileOutput << "                     VITAL STATISTICS" << '\n';
    m_FileOutput << "                     ================" << '\n';
	m_FileOutput << "      NUMBER OF DEGREES OF FREEDOM : " << m_nElements << '\n';
    m_FileOutput << "         STIFFNESS HALF-BAND WIDTH : " << m_nNodes << '\n';
	
	m_FileOutput << "                        FINAL RESULTS" << '\n';
	m_FileOutput << "                        =============" << '\n';

	//print the nodal displacements
	m_FileOutput << "\n\n";
	m_FileOutput << "NODAL DISPLACEMENTS\n";
	m_FileOutput << "Node"<<"       "<< "X DISP      "<<"       "<< "Y DISP      "<< "\n";
	m_FileOutput << "-------------------------------------------------------------------------------------\n";
	m_FileOutput << "    "<<"       "<< "   m        "<<"       "<< "   m        "<< "\n";
	m_FileOutput << "-------------------------------------------------------------------------------------\n";
	for (i=1; i <= m_nNodes; i++)
    {
        m_NodalResponseData(i).GetValues (fVResp);
        m_FileOutput <<setiosflags(std::ios::left)<< std::setw(4) << i << "   ";
        for (j=1; j <= DOFPN; j++)
            m_FileOutput << std::setw(14) << fVResp(j) << "  ";
        m_FileOutput << '\n';
    }

	
	//print the nodal reactions
	m_FileOutput << "\n\n";
	m_FileOutput << "NODAL REACTIONS                                        \n";
	m_FileOutput << "NODE"<<"       "<< "X REAC      "<<"       "<< "Y REAC      "<< "\n";
	m_FileOutput << "-------------------------------------------------------------------------------------\n";
	m_FileOutput << "    "<<"       "<< "   N        "<<"       "<< "   N        "<< "\n";
	m_FileOutput << "-------------------------------------------------------------------------------------\n";
	 for(i=1;i<=m_nNodes;i++)
	 {
		 m_NodalData(i).GetFixity (VFC);
		 if(VFC(1)==0||VFC(2)==0)
		 {
		 m_FileOutput << setiosflags(std::ios::left) << std::setw(4) << i << "       ";
		 if((VFC(1)==0)&&(VFC(2)==0))
				m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) << m_SF(2*i-1,1)<<  std::setw(10) << m_SF(2*i,1);
		 else if((VFC(1)==1)&&(VFC(2)==0))
			 m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) << "             "<<  std::setw(10) << m_SF(2*i,1);
		 else if((VFC(1)==0)&&(VFC(2)==1))
			 m_FileOutput <<setiosflags(std::ios::left)<<  std::setw(10) << m_SF(2*i-1,1)<<  std::setw(10) << "             ";
		 }
		 m_FileOutput <<"\n";
	 }
	
	//print the element stresses
	m_FileOutput << "\n\n";
	m_FileOutput << "ELEMENT STRESSES ( N / M ^2) MATERIAL GROUP : 1\n";
	m_FileOutput << "ELEM"<<"   "<< " Stress-X"<<"   "<< "Stress-Y"<< "   "<< "Shear Stress\n";
	m_FileOutput << "---------------------------------------------------------------------------------------------------------\n";
	  for (i=1; i <= m_nElements; i++)
    {
        m_ElementResponseData(i).GetStress (fVStress);
        m_FileOutput << std::setw(4) << i << "      ";
        for (j=1; j <= NUMSCOMP; j++) 
            m_FileOutput << std::setw(5) << fVStress(j) << "  ";
		m_FileOutput<<"\n";
	  }
	 
	//print the element strains
	m_FileOutput << "\n\n";
	m_FileOutput << "---------------------------------------------------------------------------------------------------------\n";
	m_FileOutput << "ELEMENT STRAINS FOR MATERIAL GROUP :   1                                   \n";
	m_FileOutput << "ELEM"<<"   "<< "Strain-X    "<<"   "<< "Strain-Y"<< "   "<< "Shear Strain\n";
	m_FileOutput << "---------------------------------------------------------------------------------------------------------\n";
	  for (i=1; i <= m_nElements; i++)
    {
        m_ElementResponseData(i).GetStrain (fVStrain);
        m_FileOutput << std::setw(4) << i << "      ";
        for (j=1; j <= NUMSCOMP; j++) 
            m_FileOutput << std::setw(5) << fVStrain(j) << "  ";
		m_FileOutput<<"\n";
	  }
	 
	//print the stress equivalents
	m_FileOutput << "\n\n";
	m_FileOutput << "------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	m_FileOutput << "STRESS EQUIVALENTS FOR MATERIAL GROUP : 1                                                       \n";
	m_FileOutput << "ELEM"<<"   "<< "    S11   "<<"   "<< "    S22   "<< "   "<< "  TRESCA  "<< "   "<< " VONMISES \n";
	m_FileOutput << "------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	  for (i=1; i <= m_nElements; i++)
    {
        m_ElementResponseData(i).GetOtherStress (fS11,fS22,fVonMis,fTersac);
        m_FileOutput << std::setw(4) << i << "      ";
            m_FileOutput << std::setw(5) << fS11 << "  "<<fS22<<"  "<<fTersac<<"  "<<fVonMis;
        m_FileOutput << '\n';
	}
}

void CFEAT::TerminateProgram ()
// ==================================================================
// Function: Closes the input and output files
//    Input: None
//   Output: None
// ==================================================================
{
    // close the input and output files
    m_FileInput.close ();
    m_FileOutput.close ();
}

void CFEAT::ErrorHandler (ErrorStates nCode)
// ==================================================================
// Function: Prints the error message and shuts down the program
//    Input: Error code
//   Output: None
// ==================================================================
{
    if (nCode == NUMNODES) 
        std::cerr << "Number of nodes must be >= 2.";
    else if (nCode == NUMELEMENTS) 
        std::cerr << "Number of elements must be >= 1.";
    else if (nCode == DEBUGCODE) 
        std::cerr << "Debug level must be 0 or 1.";
    else if (nCode == NODENUMBER) 
        std::cerr << "Invalid node number";
    else if (nCode == ELEMENTNUMBER) 
        std::cerr << "Invalid element number";
    else if (nCode == EPROPERTY) 
        std::cerr << "Element property must be positive.";
    else if (nCode == YOUNGSMODULUS) 
        std::cerr << "Modulus of elasticity must be positive.";
    else if (nCode == NODALFIXITY) 
        std::cerr << "Invalid fixity condition.";
    else if (nCode == UNSTABLE) 
        std::cerr << "Unstable structure.";
    else if (nCode == INPUTFILE) 
        std::cerr << "Unable to open specified input file.";
    else if (nCode == OUTPUTFILE) 
        std::cerr << "Unable to open specified output file.";
    else if (nCode == INVALIDINPUT) 
        std::cerr << "Invalid input value.";
    else if (nCode == INVERSEERROR)
		std::cerr << "Determinant is Zero.";

    else if (nCode == COMMANDLINE) 
    {
        std::cerr << "Invalid number of command line arguments.";
        std::cerr << "Type: feat inputfilename outputfilename\n";
    }
    if (nCode >= NUMNODES && nCode <= NODALFIXITY)
        std::cerr << '\n' << "Error in input file line : "
                  << m_nLineNumber;

    std::cerr << '\n';
    
    TerminateProgram();
    exit (1);
}