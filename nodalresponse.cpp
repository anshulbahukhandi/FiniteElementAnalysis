/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#include "FEAT.h"

CNodalResponse::CNodalResponse ()
// ==================================================================
// Function: Constructor
//    Input: None
//   Output: None
// ==================================================================
{
    m_fVResp.SetSize (DOFPN); m_fVResp.Set (0.0f);
}

CNodalResponse::~CNodalResponse ()
// ==================================================================
// Function: Destructor
//    Input: None
//   Output: None
// ==================================================================
{
}

void CNodalResponse::GetValues (CVector<float>& fVResp) const
// ==================================================================
// Function: Gets the nodal response
//    Input: None
//   Output: Vector of element response
// ==================================================================
{
    fVResp = m_fVResp;
}

void CNodalResponse::SetValues (const CVector<float>& fVResp)
// ==================================================================
// Function: Sets the nodal response
//    Input: Vector of element response
//   Output: None
// ==================================================================
{
    m_fVResp = fVResp;
}

