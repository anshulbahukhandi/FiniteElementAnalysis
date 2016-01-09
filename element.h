/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_ELEMENT_H__
#define __RAJAN_ELEMENT_H__

#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"

class CElement
{
    public:
        CElement ();
        ~CElement ();

		enum type {PLANESTRESS,PLANESTRAIN};
        // accessor functions
        void GetNodes (CVector<int>& nVList) const;
        void GetMaterialProperty (CVector<float>& fV) const;
		void Getstressstrain(type & nVNBC) const;

        // modifier functions
        void SetNodes (const CVector<int>& nVList);
		void SetMaterialProperty (const CVector<float>& fV);
		void Setstressstrain(const type & nVNBC);

    private:
        CVector<int>   m_nENodes;  // vector of element nodes
        CVector<float> m_fVMatP;   // vector of material properties
		CVector<float> m_fMatGp;
		type m_nplane;
};

#endif