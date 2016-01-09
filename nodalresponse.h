/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_NODALRESPONSE_H__
#define __RAJAN_NODALRESPONSE_H__

#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"

class CNodalResponse
{
    public:
        CNodalResponse ();
        ~CNodalResponse ();

        // accessor functions
        void GetValues (CVector<float>&) const;

        // modifier functions
        void SetValues (const CVector<float>&);

    private:
        CVector<float> m_fVResp; // vector of nodal response
};

#endif