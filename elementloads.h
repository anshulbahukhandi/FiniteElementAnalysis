
#ifndef __ELEMENTLOADS_H__
#define __ELEMENTLOADS_H__

#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"
#include "constants.h"

class CElementLoads
{
    public:
        CElementLoads ();                         // default ctor
        CElementLoads (const CElementLoads& ELO); // copy ctor
        ~CElementLoads ();                        // dtor

        // accessor functions
		void GetElementalload ( int & E,float & fv1,float & fv2) const;

        // modifier functions
		void SetElementalload (const int & E, const float & fv1, const float & fv2);
    private:
        int    m_nE;	// element number
        float  m_fV1;	// distance from start node
                            // or load intensity at start node
        float  m_fV2;	// load value or
                            // load intensity at end node
		
};		

#endif