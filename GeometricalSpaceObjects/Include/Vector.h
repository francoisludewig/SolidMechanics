#pragma once

#include <iostream>
#include <iomanip>
#include "Formatter/VectorFormatter.h"
#include "Parser/VectorParser.h"

#include "/usr/local/include/gmp.h"
#include "mpreal.h"

#define type T

namespace GeometricalSpaceObjects {
	
	template <class T>
	class Vector final {
	public:
		Vector(){
			componantX = componantY = componantZ = 0.0;
		}
		
		Vector(const T &x,const T &y, const T &z):componantX(x), componantY(y), componantZ(z){
		}
		
		virtual ~Vector() {}
		
		static Vector NullVector(){ static Vector nullVector(0,0,0); return nullVector; }
		
		void SetComponants(const T &x, const T &y, const T &z){
			this->componantX = x;
			this->componantY = y;
			this->componantZ = z;
		}
		
		T Norme() const{
			return sqrt(componantX*componantX+componantY*componantY+componantZ*componantZ);
		}
		
		void Normalize(){
			T n = Norme();
			if(n != 0){
				componantX /= n;
				componantY /= n;
				componantZ /= n;
			}
		}
		
		T ScalarProduct(const Vector & b) const{
			return (this->componantX*b.ComponantX()+this->componantY*b.ComponantY()+this->componantZ*b.ComponantZ());
		}
		
		Vector CrossProduct(const Vector & b) const{
			return std::move(Vector(componantY*b.componantZ-componantZ*b.componantY,
															componantZ*b.componantX-componantX*b.componantZ,
															componantX*b.componantY-componantY*b.componantX));
		}
		
		
		Vector Product(const T & b) const{
			return std::move(Vector(this->componantX*b, this->componantY*b, this->componantZ*b));
		}
		
		Vector Division(const T & b) const{
			return std::move(Vector(this->componantX/b, this->componantY/b, this->componantZ/b));
		}
		
		Vector Sum(const Vector & b) const{
			return std::move(Vector(this->componantX + b.componantX,
															this->componantY + b.componantY,
															this->componantZ + b.componantZ));
		}
		
		Vector Difference(const Vector & b) const{
			return std::move(Vector(this->componantX - b.componantX,
															this->componantY - b.componantY,
															this->componantZ - b.componantZ));
		}
		
		
		T ComponantX() const{return componantX;};
		T ComponantY() const{return componantY;};
		T ComponantZ() const{return componantZ;};
		
		void ComponantX(T x) {this->componantX = x;};
		void ComponantY(T y) {this->componantY = y;};
		void ComponantZ(T z) {this->componantZ = z;};
		
		void operator+=(const Vector & a){
			componantX += a.componantX;
			componantY += a.componantY;
			componantZ += a.componantZ;
		}
		
		void operator-=(const Vector & a){
			componantX -= a.componantX;
			componantY -= a.componantY;
			componantZ -= a.componantZ;
		}
		
		void operator*=(const T & a){
			componantX *= a;
			componantY *= a;
			componantZ *= a;
		}
		
		void operator/=(const T & a){
			componantX /= a;
			componantY /= a;
			componantZ /= a;
		}
		
		
		T operator*(const Vector &b) const{
			return this->ScalarProduct(b);
		}
		
		Vector operator^(const Vector &b) const{
			return this->CrossProduct(b);
		}
		
		Vector operator*(const T &b) const{
			return this->Product(b);
		}
		
		Vector operator/(const T &b) const{
			return this->Division(b);
		}
		
		Vector operator+(const Vector &b) const{
			return this->Sum(b);
		}
		
		Vector operator-(const Vector &b) const{
			return this->Difference(b);
		}
		
		std::string Format(VectorFormatter<T>* vectorFormatter) {vectorFormatter->Format(*this);}
		
		void Parse(VectorParser<T>* vectorParser, const std::string& str) {*this = vectorParser->Parse(str);}
		
	private:
		T componantX,componantY,componantZ;
	};
	
	template <typename T>
	Vector<T> operator*(const T & b, const Vector<T> & a){
		return a.Product(b);
	}
	
	template<class T>
	bool operator== (const Vector<T> &vector1, const Vector<T> &vector2) {
		if(fabs(vector1.ComponantX()-vector2.ComponantX()) > std::numeric_limits<T>::epsilon())
			return false;
		
		if(fabs(vector1.ComponantY()-vector2.ComponantY()) > std::numeric_limits<T>::epsilon())
			return false;
		
		if(fabs(vector1.ComponantZ()-vector2.ComponantZ()) > std::numeric_limits<T>::epsilon())
			return false;
		
		return true;
	}
	
	template<class T>
	bool operator!= (const Vector<T> &vector1, const Vector<T> &vector2){
		return !(vector1==vector2);
	}
	
}

template<class T> inline std::ostream & operator << (std::ostream & out, const GeometricalSpaceObjects::Vector<T> & a) {
	out.precision(16);
	out << std::scientific << a.ComponantX() << "\t" << a.ComponantY() << "\t" << a.ComponantZ();
	return out;
}

template<class T> inline std::istream & operator >> (std::istream & in, GeometricalSpaceObjects::Vector<T> & a){
	T x,y,z;
	in >> x >> y >> z;
	a.SetComponants(x,y,z);
	return in;
}


