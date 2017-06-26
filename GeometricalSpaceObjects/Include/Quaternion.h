#pragma once

#include <iostream>
#include <iomanip>
#include "Formatter/QuaternionFormatter.h"
#include "Parser/QuaternionParser.h"

#include "/usr/local/include/gmp.h"
#include "mpreal.h"


namespace GeometricalSpaceObjects {
	
	template<class T>
	class Vector;
	
	template<class T>
	class Quaternion final {
	public:
		Quaternion(){
			componantReal = 1.0;
			componantI = componantJ = componantK = 0.0;
		}

		Quaternion(const T & q0, const T & q1, const T & q2, const T & q3){
			this->componantReal = q0;
			this->componantI = q1;
			this->componantJ = q2;
			this->componantK = q3;
		}

		Quaternion(const Vector<T> & w){
			T a = w.Norme();
			if(a != 0){
				T sa = sin(a/2);
				T ca = cos(a/2);
				componantReal = ca;
				componantI = w.ComponantX()/a*sa;
				componantJ = w.ComponantY()/a*sa;
				componantK = w.ComponantZ()/a*sa;
			}
		}

		~Quaternion() {}
		
		T ComponantReal() const{return componantReal;}
		T ComponantI() const{return componantI;}
		T ComponantJ() const{return componantJ;}
		T ComponantK() const{return componantK;}
		T Norme(){
			return sqrt(componantReal*componantReal+componantI*componantI+componantJ*componantJ+componantK*componantK);
		}

		void Normalize(){
			T n = Norme();
			if(n != 0){
				componantReal /= n;
				componantI /= n;
				componantJ /= n;
				componantK /= n;
			}
		}

		void SetComponants(const T & q0, const T & q1, const T & q2, const T & q3){
			this->componantReal = q0;
			this->componantI = q1;
			this->componantJ = q2;
			this->componantK = q3;
		}

		
		Quaternion Product(const Quaternion & b) const{
			return std::move(Quaternion(this->componantReal*b.ComponantReal() - this->componantI*b.ComponantI() - this->componantJ*b.
																	ComponantJ() - this->componantK*b.ComponantK(),
																	this->componantReal*b.ComponantI() + this->componantI*b.ComponantReal() - this->componantJ*b.
																	ComponantK() + this->componantK*b.ComponantJ(),
																	this->componantReal*b.ComponantJ() + this->componantI*b.ComponantK() + this->componantJ*b.
																	ComponantReal() - this->componantK*b.ComponantI(),
																	this->componantReal*b.ComponantK() - this->componantI*b.ComponantJ() + this->componantJ*b.
																	ComponantI() + this->componantK*b.ComponantReal()));
		}

		Quaternion Sum(const Quaternion & b) const{
			return std::move(Quaternion(this->componantReal + b.componantReal,
																	this->componantI + b.componantI,
																	this->componantJ + b.componantJ,
																	this->componantK + b.componantK));
		}

		Quaternion Diff(const Quaternion & b) const{
			return std::move(Quaternion(this->componantReal - b.componantReal,
																	this->componantI - b.componantI,
																	this->componantJ - b.componantJ,
																	this->componantK - b.componantK));
		}
		Quaternion operator*(const Quaternion &b) const{
			return this->Product(b);
		}

		Quaternion operator+(const Quaternion &b) const{
			return this->Sum(b);
		}

		Quaternion operator-(const Quaternion &b) const{
			return this->Diff(b);
		}

		Quaternion operator~(){
			return std::move(Quaternion(this->componantReal,
																	-this->componantI,
																	-this->componantJ,
																	-this->componantK));
		}

		void operator*=(const Quaternion & b){
			T p0 = this->componantReal*b.ComponantReal() - this->componantI*b.ComponantI() - this->componantJ*b.ComponantJ() - this->componantK*b.ComponantK();
			T p1 = this->componantReal*b.ComponantI() + this->componantI*b.ComponantReal() - this->componantJ*b.ComponantK() + this->componantK*b.ComponantJ();
			T p2 = this->componantReal*b.ComponantJ() + this->componantI*b.ComponantK() + this->componantJ*b.ComponantReal() - this->componantK*b.ComponantI();
			T p3 = this->componantReal*b.ComponantK() - this->componantI*b.ComponantJ() + this->componantJ*b.ComponantI() + this->componantK*b.ComponantReal();
			componantReal = p0;
			componantI = p1;
			componantJ = p2;
			componantK = p3;
		}

		void operator+=(const Quaternion & a){
			componantReal += a.ComponantReal();
			componantI += a.ComponantI();
			componantJ += a.ComponantJ();
			componantK += a.ComponantK();
		}

		void operator-=(const Quaternion & a){
			componantReal -= a.ComponantReal();
			componantI -= a.ComponantI();
			componantJ -= a.ComponantJ();
			componantK -= a.ComponantK();
		}

		std::string Format(QuaternionFormatter<T>* quaternionFormatter) {quaternionFormatter->Format(*this);}
		
		void Parse(QuaternionParser<T>* quaternionParser, const std::string& str) {*this = quaternionParser->Parse(str);}
		
		
	private:
		T componantReal,componantI,componantJ,componantK;
	};
	
}

template<class T> inline std::ostream & operator << (std::ostream & out, const GeometricalSpaceObjects::Quaternion<T> & a){
	out.precision(16);
	out << std::scientific << a.ComponantReal() << "\t" << a.ComponantI() << "\t" << a.ComponantJ() << "\t" << a.ComponantK();
	return out;
}

template<class T> inline std::istream & operator >> (std::istream & in, GeometricalSpaceObjects::Quaternion<T> & a){
	T q0,q1,q2,q3;
	in >> q0 >> q1 >> q2 >> q3;
	a.SetComponants(q0,q1,q2,q3);
	return in;
}
