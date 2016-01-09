/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#include "FEAT.h"

CElementResponse::CElementResponse ()
// ==================================================================
// Function: Constructor
//    Input: None
//   Output: None
// ==================================================================
{
    m_fVStrain.SetSize (NUMSCOMP);
    m_fVStress.SetSize (NUMSCOMP);
}

CElementResponse::~CElementResponse ()
// ==================================================================
// Function: Destructor
//    Input: None
//   Output: None
// ==================================================================
{
}

void CElementResponse::GetStrain (CVector<float>& fVS) const
// ==================================================================
// Function: Gets the element strains
//    Input: None
//   Output: Vector of element strains
// ==================================================================
{
    fVS = m_fVStrain;
}
void CElementResponse::GetStress (CVector<float>& fVS) const
// ==================================================================
// Function: Gets the element strains
//    Input: None
//   Output: Vector of element strains
// ==================================================================
{
    fVS = m_fVStress;
}

void CElementResponse::GetOtherStress (float& fS11,float &fS22,float &fVonMis,float &fTresca) const
// ==================================================================
// Function: Gets the element stresses
//    Input: None
//   Output: Vector of element stresses
// ==================================================================
{
    fS11 = m_fS11;
	fS22 = m_fS22;
	fVonMis = m_fVonMis;
	fTresca = m_fTresca;
}

void CElementResponse::SetStrain (const CVector<float>& fVS)
// ==================================================================
// Function: Sets the element strains
//    Input: Vector of element strains
//   Output: None
// ==================================================================
{
    m_fVStrain = fVS;
}
void CElementResponse::SetStress (const CVector<float>& fVS)
// ==================================================================
// Function: Sets the element strains
//    Input: Vector of element strains
//   Output: None
// ==================================================================
{
    m_fVStress = fVS;
}
void CElementResponse::SetOtherStress (float &fS11,float &fS22,float &fVonMis,float &fTresca)
// ==================================================================
// Function: Sets the element stresses
//    Input: None
//   Output: Vector of element stresses
// ==================================================================
{
    m_fS11 = fS11;
	m_fS22 = fS22;
	m_fVonMis = fVonMis;
	m_fTresca = fTresca;
}

