/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#include "FEAT.h"

CElement::CElement ()
// ==================================================================
// Function: Constructor
//    Input: None
//   Output: None
// ==================================================================
{
    m_nENodes.SetName ("m_nENodes");
    m_nENodes.SetSize (NUMENODES); m_nENodes.Set(0);
	m_fMatGp.SetSize(NUMEP);
    m_fVMatP.SetSize (NUMMATP);    m_fVMatP.Set (0.0f);

}

CElement::~CElement ()
// ==================================================================
// Function: Destructor
//    Input: None
//   Output: None
// ==================================================================
{
}

void CElement::GetNodes (CVector<int>& nVList) const
// ==================================================================
// Function: Gets the list of element nodes
//    Input: None
//   Output: Vector of element nodes
// ==================================================================
{
    for (int i=1; i <= NUMENODES; i++)
        nVList(i) = m_nENodes(i);
}


void CElement::SetNodes (const CVector<int>& nVList)
// ==================================================================
// Function: Sets the list of element nodes
//    Input: Vector of element nodes
//   Output: None
// ==================================================================
{
    for (int i=1; i <= NUMENODES; i++)
        m_nENodes(i) = nVList(i);
}

void CElement::SetMaterialProperty(const CVector<float>& nVList)
{
	for (int i=1; i <= NUMEP; i++)
         m_fMatGp(i) = nVList(i);
}

void CElement::GetMaterialProperty(CVector<float>& nVList) const
{
	for (int i=1; i <= NUMEP; i++)
       nVList(i) = m_fMatGp(i);
}

void CElement::Setstressstrain(const type & nVNBC)
{
	m_nplane = nVNBC;
	
}

void CElement::Getstressstrain(type & nVNBC) const
{
	nVNBC=m_nplane;
}