/***********************************************
A  FEA  Template  Program
Copyright(c) 2004-11, S. D. Rajan
All rights reserved
***********************************************/
#include <sstream>
#include "FEAT.h"
#include "F:\ASU Courses\cee532\Library\Library\MatToolBox.h"
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\NumericalIntegration.h"

CMatrix<double> DB(3,8),DEO(3,1),DD(3,3);


int CFEAT::k1DC0LElement (int nE, CVector<int>& nVNList,
						  CMatrix<double>& dMk,CMatrix<double> & dTF)
// ==================================================================
// Function: Generates k for the 1D-C0 linear element
//    Input: Element number
//   Output: Element stiffness matrix in dMk
//           a return value of 0 means no error
// ==================================================================
{
          
    CVector<float> fVMatP(NUMEP);     // young's modulus
    CVector<float> fVCoords1(PDIM);     // nodal coordinates (node 1)
    CVector<float> fVCoords2(PDIM);     // nodal coordinates (node 2)
	 CVector<float> fVCoords3(PDIM);    // nodal coordinates (node 3)
	 CVector<float> fVCoords4(PDIM);    // nodal coordinates (node 4)
	 CVector<int> nVN(NUMENODES);
	 CElement::type VC;
    // get element data
    m_ElementData(nE).GetMaterialProperty(fVMatP);
    m_ElementData(nE).GetNodes(nVN);
	m_ElementData(nE).Getstressstrain(VC);
    m_NodalData(nVN(1)).GetCoords (fVCoords1);
    m_NodalData(nVN(2)).GetCoords (fVCoords2);
	m_NodalData(nVN(3)).GetCoords (fVCoords3);
	m_NodalData(nVN(4)).GetCoords (fVCoords4);
	float x = 0.0f;
    CGaussQuad GQ;          // for numerical integration via G-Q
    const int NORDER = 2;   // numerical integration order
    const int NGP    = 2;   // # of G-Q points
    CVector<double> PHI(NUMENODES);    // shape functions
    CVector<double> DPHDXI(NUMENODES); // derivatives of shape functions
    CVector<double> DPHDX(NUMENODES);  // derivatives of shape functions (wrt x)
	CMatrix<double> DJAC(NORDER,NORDER);
	CMatrix<double> DJACIN(NORDER,NORDER);
	CMatrix<double> DO(3,4),DN(4,8),DBT(8,3),DT(8,3),DK(8,8),DTE(8,1),DTI(3,3),DTT(8,3);
	double delta=0.0,DETJAC;
	float T;
	DB.Set(0.0);
	DN.Set(0.0);
	DO.Set(0.0);
	DBT.Set(0.0);
	DT.Set(0.0);
	DEO.Set(0.0);
	DTE.Set(0.0);
	DD.Set(0.0);
	DK.Set(0.0);
	DTI.Set(0.0);
	DTT.Set(0.0);
    int i, j;

    for (int NG=1; NG <= NGP; NG++)
    {
		for (int NN=1; NN<=NGP; NN++)
		{
        // gauss-point and weight
        double PXI = GQ.GetLocation (NORDER, NG);
		double PNI = GQ.GetLocation (NORDER, NN);
        double WG1  = GQ.GetWeight (NORDER, NG);
		double WG2  = GQ.GetWeight (NORDER, NN);
        // shape functions
        PHI(1)=0.25*(1.0-PXI)*(1.0-PNI);
		PHI(2)=0.25*(1.0+PXI)*(1.0-PNI);
		PHI(3)=0.25*(1.0+PXI)*(1.0+PNI);
		PHI(4)=0.25*(1.0-PXI)*(1.0+PNI);

		double a=fVCoords1(1)-fVCoords2(1)+fVCoords3(1)-fVCoords4(1);
		double a1=-fVCoords1(1)+fVCoords2(1)+fVCoords3(1)-fVCoords4(1);
		double a2=-fVCoords1(1)-fVCoords2(1)+fVCoords3(1)+fVCoords4(1);
		double b=fVCoords1(2)-fVCoords2(2)+fVCoords3(2)-fVCoords4(2);
		double b1=-fVCoords1(2)+fVCoords2(2)+fVCoords3(2)-fVCoords4(2);
		double b2=-fVCoords1(2)-fVCoords2(2)+fVCoords3(2)+fVCoords4(2);
		
        // compute jacobian
		DJAC(1,1)=(PNI*a+a1)*0.25;
		DJAC(1,2)=(PNI*b+b1)*0.25;
		DJAC(2,1)=(PXI*a+a2)*0.25;
		DJAC(2,2)=(PXI*b+b2)*0.25;
		
		DETJAC=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1);
		
		/*if(DETJAC==0)
			ErrorHandler(INVERSEERROR);*/
        // compute inverse of jacobian
		DJACIN(1,1)=DJAC(2,2)/DETJAC;
		DJACIN(2,2)=DJAC(1,1)/DETJAC;
		if(DJAC(1,2)!=0)
		DJACIN(1,2)=-DJAC(1,2)/DETJAC;
		else
		DJACIN(1,2)=DJAC(1,2)/DETJAC;
		if(DJAC(2,1)!=0)
		DJACIN(2,1)=-DJAC(2,1)/DETJAC;
		else
		DJACIN(2,1)=DJAC(2,1)/DETJAC;
		DO(1,1)=DJACIN(1,1);
		DO(1,2)=DJACIN(1,2);
		DO(2,3)=DJACIN(2,1);
		DO(2,4)=DJACIN(2,2);
		DO(3,1)=DJACIN(2,1);
		DO(3,2)=DJACIN(2,2);
		DO(3,3)=DJACIN(1,1);
		DO(3,4)=DJACIN(1,2);

		DN(1,1)=DN(3,2)= -(1-PNI);
		DN(1,3)=DN(3,4)= -DN(1,1);
		DN(1,5)=DN(3,6)=  (1+PNI);
		DN(1,7)=DN(3,8)= -DN(1,5);
		DN(2,1)=DN(4,2)= -(1-PXI);
		DN(2,3)=DN(4,4)= -(1+PXI);
		DN(2,5)=DN(4,6)= -DN(2,3);
		DN(2,7)=DN(4,8)= -DN(2,1);
		if (VC==0)
		{
		DD(1,1)=DD(2,2)=fVMatP(2)/(1-fVMatP(3)*fVMatP(3));
		DD(1,2)=DD(2,1)=fVMatP(3)*fVMatP(2)/(1-fVMatP(3)*fVMatP(3));
		DD(3,3)=0.5*(1-fVMatP(3))*fVMatP(2)/(1-fVMatP(3)*fVMatP(3));
		}
		else
		{
		DD(1,1)=DD(2,2)=(1-fVMatP(3))*fVMatP(2)/((1+fVMatP(3))*(1-2*fVMatP(3)));
		DD(1,2)=DD(2,1)=fVMatP(3)*fVMatP(2)/((1+fVMatP(3))*(1-2*fVMatP(3)));
		DD(3,3)=(0.5-fVMatP(3))*fVMatP(2)/((1+fVMatP(3))*(1-2*fVMatP(3)));
		}
		Multiply(DO,DN,DB);
		for(int h=1;h<=4;h++)
			{
				m_NodalData(nVN(h)).GetTemp(T);
				delta+=PHI(h)*T;
			}
			if(VC==0)
				DEO(1,1)=DEO(2,1)=fVMatP(4)*delta;
			else
				DEO(1,1)=DEO(2,1)=(1+fVMatP(3))*fVMatP(4)*delta;
			
        // compute stiffness at gauss point
        for (i=1; i <= 3; i++)
			for (j=1; j <= 3; j++)
				DTI(i,j) += WG1*WG2*DD(i,j)*fVMatP(1)*DETJAC;	
		Transpose(DB,DBT);
		Multiply(DBT,DTI,DT);
		Multiply(DBT,DD,DTT);
		Multiply(DTT,DEO,DTE);
		Multiply(DT,DB,DK);
		for(int k=1;k<=2*NUMENODES;k++)
		{
			dTF(k,1) += WG1*WG2*DTE(i,1)*fVMatP(1);  
			for(int l=1;l<=2*NUMENODES;l++)
				dMk(k,l)+=DK(k,l);
		}
	}
}


	
    // debug?
    return 0;
}

void CFEAT::s1DC0LElement (int nE) 
// ==================================================================
// Function: Computes the element strain and stress
//           for the 1D-C0 linear element
//    Input: Element number
//   Output: None
// ==================================================================
{
    CVector<int>   nVNL(NUMENODES);  // list of element nodes
  
    CVector<float> fVMatP(NUMEP);     // young's modulus
    CVector<float> fVCoords1(PDIM);     // nodal coordinates (node 1)
    CVector<float> fVCoords2(PDIM);     // nodal coordinates (node 2)
	CVector<float> fVCoords3(PDIM);    // nodal coordinates (node 3)
	CVector<float> fVCoords4(PDIM);    // nodal coordinates (node 4)
    CVector<float> fVDisp1(DOFPN);      // nodal displacements (node 1)
    CVector<float> fVDisp2(DOFPN);      // nodal displacements (node 2)
    CVector<float> fVStrain(NUMSCOMP);  // element strain
    CVector<float> fVStress(NUMSCOMP);  // element stress
	CMatrix<double> DE(3,1),DS(3,1),DSR(3,1),DDI(8,1);
	float fVonMis,fTresca,fS11,fS22;
	DE.Set(0.0);
	DS.Set(0.0);
	DSR.Set(0.0);
    // get element data

    m_ElementData(nE).GetMaterialProperty (fVMatP);
    m_ElementData(nE).GetNodes (nVNL);
    m_NodalData(nVNL(1)).GetCoords (fVCoords1);
    m_NodalData(nVNL(2)).GetCoords (fVCoords2);
	m_NodalData(nVNL(3)).GetCoords (fVCoords3);
    m_NodalData(nVNL(4)).GetCoords (fVCoords4);
    m_NodalResponseData(nVNL(1)).GetValues (fVDisp1);
    m_NodalResponseData(nVNL(2)).GetValues (fVDisp2);
	int l=1;
	for(int i=1;i<=NUMENODES;i++)
	{
		DDI(l,1)=m_SND(2*nVNL(i)-1,1);
		l++;
		DDI(l,1)=m_SND(2*nVNL(i),1);
		l++;
	}
	Multiply(DB,DDI,DE);
	Subtract(DE,DEO,DS);
   Multiply(DD,DS,DSR);
    // strain
   for(int i=1;i<=3;i++)
    fVStrain(i) = DE(i,1);

    // stress
    for(int i=1;i<=3;i++)
    fVStress(i) = DSR(i,1);
	fS11=((DSR(1,1)+DSR(2,1))/2)+0.5*sqrt(pow(DSR(1,1)-DSR(2,1),2)+4*pow(DSR(3,1),2));
	fS22=((DSR(1,1)+DSR(2,1))/2)-0.5*sqrt(pow(DSR(1,1)-DSR(2,1),2)+4*pow(DSR(3,1),2));
	fVonMis=sqrt(pow(fS11,2)+pow(fS22,2)-fS11*fS22);
	fTresca=abs(fS11-fS22)/2;
    // update model with the computed values
    m_ElementResponseData(nE).SetStrain (fVStrain);
	m_ElementResponseData(nE).SetStress(fVStress);
	m_ElementResponseData(nE).SetOtherStress (fS11,fS22,fVonMis,fTresca);
}