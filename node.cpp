/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#include "FEAT.h"

CNode::CNode ()
// ==================================================================
// Function: Constructor
//    Input: None
//   Output: None
// ==================================================================
{
    m_fVCoor.SetSize(PDIM);    m_fVCoor.Set (0.0f);
	m_nVNBC.SetSize(DOFPN);    m_nVNBC.Set(CNode::FREE);
    m_fVForces.SetSize(DOFPN); m_fVForces.Set(0.0f); 
	m_Temp=0.0f; m_fDisp.SetSize(DOFPN);m_fDisp.Set(0.0f);
}

CNode::~CNode ()
// ==================================================================
// Function: Destructor
//    Input: None
//   Output: None
// ==================================================================
{
}

void CNode::GetCoords (CVector<float>& fVCoor) const
// ==================================================================
// Function: Gets the nodal coordinates
//    Input: None
//   Output: Vector of nodal coordinates
// ==================================================================
{
	for(int i=1;i<=2;i++)
    fVCoor(i) = m_fVCoor(i);
}

void CNode::SetCoords (const CVector<float>& fVCoor)
// ==================================================================
// Function: Sets the nodal coordinates
//    Input: Vector of nodal coordinates
//   Output: None
// ==================================================================
{
	for(int i=1;i<=2;i++)
    m_fVCoor(i) = fVCoor(i);
}

void CNode::GetFixity (CVector<NBConditions>& nVNBC) const
// ==================================================================
// Function: Gets the nodal fixity conditions
//    Input: None
//   Output: Vector of fixity conditions
// ==================================================================
{
	for(int i=1;i<=2;i++)
    nVNBC(i) = m_nVNBC(i);
}

void CNode::SetFixity (const CVector<NBConditions>& nVNBC)
// ==================================================================
// Function: Sets the nodal fixity conditions
//    Input: Vector of fixity conditions
//   Output: None
// ==================================================================
{
	for(int i=1;i<=2;i++)
    m_nVNBC(i) = nVNBC(i);
}

void CNode::GetLoads (CVector<float>& fVForces) const
// ==================================================================
// Function: Gets the nodal forces (applied loads)
//    Input: None
//   Output: Vector of nodal forces
// ==================================================================
{
    fVForces = m_fVForces;
}

void CNode::SetLoads (const CVector<float>& fVForces)
// ==================================================================
// Function: Sets the nodal forces (applied loads)
//    Input: Vector of nodal forces
//   Output: None
// ==================================================================
{
    m_fVForces = fVForces;
}

void CNode::GetDisps (CVector<float>& fVDisp) const
// ==================================================================
// Function: Gets the nodal coordinates
//    Input: None
//   Output: Vector of nodal coordinates
// ==================================================================
{
	for(int i=1;i<=2;i++)
    fVDisp(i) = m_fDisp(i);
}

void CNode::SetDisps (const CVector<float>& fVDisp) 
// ==================================================================
// Function: Gets the nodal coordinates
//    Input: None
//   Output: Vector of nodal coordinates
// ==================================================================
{
	for(int i=1;i<=2;i++)
    m_fDisp(i) = fVDisp(i);
}

void CNode::GetTemp (float & fVForces) const
// ==================================================================
// Function: Gets the nodal forces (applied loads)
//    Input: None
//   Output: Vector of nodal forces
// ==================================================================
{
    fVForces = m_Temp  ;
}

void CNode::SetTemp (const float & fVForces)
// ==================================================================
// Function: Sets the nodal forces (applied loads)
//    Input: Vector of nodal forces
//   Output: None
// ==================================================================
{
     m_Temp = fVForces ;
}