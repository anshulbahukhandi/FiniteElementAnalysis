/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_ELEMENTRESPONSE_H__
#define __RAJAN_ELEMENTRESPONSE_H__

#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"

class CElementResponse
{
    public:
        CElementResponse ();
        ~CElementResponse ();

        // accessor functions
        void GetStrain (CVector<float>& fVS) const;
		void GetStress (CVector<float>& fVS) const;
        void GetOtherStress (float& fS11,float &fS22,float &fVonMis,float &fTresca) const;

        // modifier functions
        void SetStrain (const CVector<float>& fVS);
		void SetStress (const CVector<float>& fVS);
        void SetOtherStress (float& fS11, float & fS22, float &fVonMis, float &fTresca);

    private:
        CVector<float> m_fVStrain;  // vector of element strains
        CVector<float> m_fVStress;  // vector of element stresses
		float m_fS11,m_fS22,m_fVonMis,m_fTresca;
};

#endif