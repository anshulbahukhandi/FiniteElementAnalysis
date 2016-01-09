
#include "ElementLoads.h"

CElementLoads::CElementLoads ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nE = 0;
    m_fV1 = m_fV2 = 0.0f;
}

CElementLoads::CElementLoads (const CElementLoads& ELO)
// ---------------------------------------------------------------------------
// Function: copy constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nE = ELO.m_nE;
    m_fV1 = ELO.m_fV1;
    m_fV2 = ELO.m_fV2;
}

CElementLoads::~CElementLoads ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementLoads::GetElementalload(int & E,float & fv1,float & fv2)const
{
	E = m_nE;
	fv1=m_fV1;
	fv2=m_fV2;
}

void CElementLoads::SetElementalload(const int & E, const float & fv1, const float & fv2)
{
	m_nE=E;
	m_fV1=fv1;
	m_fV2=fv2;
}