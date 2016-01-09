/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_NODE_H__
#define __RAJAN_NODE_H__

#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"

class CNode
{
    public:
        CNode ();
        ~CNode ();

        enum NBConditions {KNOWN, FREE};

        // accessor functions
        void GetCoords (CVector<float>&) const;
        void GetFixity (CVector<NBConditions>&) const;
		void GetDisps (CVector<float>&) const;
        void GetLoads  (CVector<float>&) const;
        void GetTemp  (float &) const;

        // modifier functions
        void SetCoords (const CVector<float>&);
		void SetDisps (const CVector<float>&);
        void SetFixity (const CVector<NBConditions>&);
        void SetLoads  (const CVector<float>&);
        void SetTemp  (const float &);


    private:
        CVector<float>        m_fVCoor;   // vector of nodal coordinates
        CVector<NBConditions> m_nVNBC;    // vector of nodal boundary conditions
        CVector<float>        m_fVForces; // vector of nodal forces
		CVector<float>		  m_fDisp;
		float m_Temp;
};

#endif  