#pragma once

#include <iostream>
#include <iomanip>
#include "Quaternion.h"
#include "Vector.h"
#include "Point.h"
#include "VectorsQuaternionConverter.h"
#include "Formatter/BasisFormatter.h"
#include "Parser/BasisParser.h"
#include "VectorsQuaternionConverter.h"

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Basis {
	public:
		Basis():axisX(1,0,0),axisY(0,1,0),axisZ(0,0,1),origin() {}
		Basis(const Point<T> & o, const Quaternion<T> & q):origin(o), orientation(q) {
			this->vQc.ConvertQuaternionIntoVectors(q,axisX,axisY,axisZ);
		}
		
		~Basis() {}
		
		Point<T> Origin() const {return this->origin;}
		Vector<T> AxisX() const {return this->axisX;}
		Vector<T> AxisY() const {return this->axisY;}
		Vector<T> AxisZ() const {return this->axisZ;}
		Quaternion<T> Orientation() const {return this->orientation;}
		
		void AxisX(const Vector<T> & e1){
			this->axisX = e1;
			this->axisX.Normalize();
			this->ConstructAxisYAndZFromX();
		}
		
		void Origin(const Point<T> & o) {this->origin = o;}
		
		void Orientation(const Quaternion<T> & q) {
			this->orientation = q;
			this->vQc.ConvertQuaternionIntoVectors(q,axisX,axisY,axisZ);
		}
		
		void Rotate(const Quaternion<T> & q) {
			this->orientation *= q;
			this->vQc.ConvertQuaternionIntoVectors(this->orientation,axisX,axisY,axisZ);
		}
		
		void Translate(const Vector<T> & o){
			this->origin += o;
		}
		
		void Local(Vector<T> & a) const{
			T x,y,z;
			x = a*axisX;
			y = a*axisY;
			z = a*axisZ;
			a.SetComponants(x,y,z);
		}
		
		void Global(Vector<T> & a) const{
			a = a.ComponantX()*axisX + a.ComponantY()*axisY + a.ComponantZ()*axisZ;
		}
		
		Point<T> Local(const Point<T> & a) const{
			Vector<T> v = a-origin;
			return Point<T>(axisX*v,axisY*v,axisZ*v);
		}
		
		Point<T> Global(const Point<T> & a) const{
			Point<T> b = origin;
			b += a.CoordinateX()*axisX + a.CoordinateY()*axisY + a.CoordinateZ()*axisZ;
			return b;
		}
		
		void LoadFromIstream(std::istream & in){
			in >> this->origin;
			in >> this->orientation;
			this->vQc.ConvertQuaternionIntoVectors(this->orientation,this->axisX,this->axisY,this->axisZ);
		}
		
		Basis<T> operator*(const Quaternion<T> & q) const{
			Basis b = *this;
			b.Rotate(q);
			return b;
		}
		
		Basis<T> operator+(const Vector<T> & o) const{
			Basis b = *this;
			b.Translate(o);
			return b;
		}
		
		Basis<T>& operator*=(const Quaternion<T> & q){
			this->Rotate(q);
			return *this;
		}
		
		Basis<T>& operator+=(const Vector<T> & o){
			this->Translate(o);
			return *this;
		}
		
		std::string Format(BasisFormatter<T>* basisFormatter) {basisFormatter->Format(*this);}
		
		void Parse(BasisParser<T>* basisParser, const std::string& str) {*this = basisParser->Parse(str);}
		
	private:
		void ConstructAxisYAndZFromX() {
			if((axisX.ComponantX() != 0 || axisX.ComponantY() != 0) || (axisX.ComponantX() != 0 || axisX.ComponantZ() != 0)  || (axisX.ComponantY() != 0 || axisX.ComponantZ() != 0)){
				axisY.ComponantX(axisX.ComponantY()*axisX.ComponantZ());
				axisY.ComponantY(axisX.ComponantX()*axisX.ComponantZ());
				axisY.ComponantZ(-2*axisX.ComponantX()*axisX.ComponantY());
			}
			else{
				if(axisX.ComponantX() == 0 || axisX.ComponantY() == 0){
					axisY.ComponantX(axisX.ComponantZ());
					axisY.ComponantY(0);
					axisY.ComponantZ(0);
				}
				else if(axisX.ComponantX() == 0 || axisX.ComponantZ() == 0){
					axisY.ComponantX(0);
					axisY.ComponantY(0);
					axisY.ComponantZ(axisX.ComponantY());
				}
				else if(axisX.ComponantY() == 0 || axisX.ComponantZ() == 0){
					axisY.ComponantX(0);
					axisY.ComponantY(axisX.ComponantX());
					axisY.ComponantZ(0);
				}
			}
			axisY.Normalize();
			axisZ = axisX^axisY;
		}
		
		Vector<T> axisX,axisY,axisZ;
		Point<T> origin;
		Quaternion<T> orientation;
		VectorsQuaternionConverter<T> vQc;
	};

}

template<class T> inline std::ostream & operator << (std::ostream & out, const GeometricalSpaceObjects::Basis<T> & a) {
	out << a.Origin() << "\n" << a.Orientation();
	return out;
}

template<class T> inline std::istream & operator >> (std::istream & in, GeometricalSpaceObjects::Basis<T> & a){
	a.LoadFromIstream(in);
	return in;
}
