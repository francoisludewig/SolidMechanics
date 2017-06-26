#pragma once

#include <iostream>
#include <iomanip>

namespace GeometricalSpaceObjects {
	template<class T>
	class Vector;
	template<class T>
	class Quaternion;
	
	template<class T>
	class VectorsQuaternionConverter {
	public:
		VectorsQuaternionConverter() {}
		~VectorsQuaternionConverter() {}

		void ConvertVectorsIntoQuaternion(const Vector<T>& e1,const Vector<T>& e2,const Vector<T>& e3, Quaternion<T>& q) const{
			T q0,q1,q2,q3;
			if(-(e2.ComponantY()+e3.ComponantZ()-e1.ComponantX()-1) < 0)
				q1 = 0.;
			else
				q1 = sqrt(-(e2.ComponantY()+e3.ComponantZ()-e1.ComponantX()-1)/4.);
			
			if(-(-e2.ComponantY()+e3.ComponantZ()+e1.ComponantX()-1) < 0)
				q2 = 0.;
			else
				q2 = sqrt(-(-e2.ComponantY()+e3.ComponantZ()+e1.ComponantX()-1)/4.);
			
			if(-(e2.ComponantY()-e3.ComponantZ()+e1.ComponantX()-1) < 0)
				q3 = 0.;
			else
				q3 = sqrt(-(e2.ComponantY()-e3.ComponantZ()+e1.ComponantX()-1)/4.);
			
			if(1.-q1*q1-q2*q2-q3*q3 >= 0)
				q0 = sqrt(1.-q1*q1-q2*q2-q3*q3);
			else
				q0 = 0.;
			
			// Cancel quadratic relative errors
			T maxq = q0*q0;
			if(q1 > maxq)maxq = q1*q1;
			if(q2 > maxq)maxq = q2*q2;
			if(q3 > maxq)maxq = q3*q3;
			T error = maxq * std::numeric_limits<T>::epsilon();
			if(q0*q0 < error)q0 = 0.;
			if(q1*q1 < error)q1 = 0.;
			if(q2*q2 < error)q2 = 0.;
			if(q3*q3 < error)q3 = 0.;
			// Sign correction
			if(e2.ComponantZ()-e3.ComponantY() < 0) q1*=-1;
			if(e3.ComponantX()-e1.ComponantZ() < 0) q2*=-1;
			if(e1.ComponantY()-e2.ComponantZ() < 0) q3*=-1;
			
			q.SetComponants(q0,q1,q2,q3);
		}

		void ConvertQuaternionIntoVectors(const Quaternion<T>& q, Vector<T>& e1, Vector<T>& e2, Vector<T>& e3) const{
			e1.SetComponants(1 - 2*q.ComponantJ()*q.ComponantJ() - 2*q.ComponantK()*q.ComponantK(),
											 2*q.ComponantI()*q.ComponantJ() + 2*q.ComponantK()*q.ComponantReal(),
											 2*q.ComponantI()*q.ComponantK() - 2*q.ComponantJ()*q.ComponantReal());
			
			e2.SetComponants(2*q.ComponantI()*q.ComponantJ() - 2*q.ComponantK()*q.ComponantReal(),
											 1 - 2*q.ComponantI()*q.ComponantI() - 2*q.ComponantK()*q.ComponantK(),
											 2*q.ComponantJ()*q.ComponantK() + 2*q.ComponantI()*q.ComponantReal());
			
			
			e3.SetComponants(2*q.ComponantI()*q.ComponantK() + 2*q.ComponantJ()*q.ComponantReal(),
											 2*q.ComponantJ()*q.ComponantK() - 2*q.ComponantI()*q.ComponantReal(),
											 1 - 2*q.ComponantJ()*q.ComponantJ() - 2*q.ComponantI()*q.ComponantI());
				}

	};
	
}
