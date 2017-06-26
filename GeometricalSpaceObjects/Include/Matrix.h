#pragma once

#include <iostream>
#include <iomanip>

#include "Quaternion.h"
#include "Vector.h"
#include "Formatter/MatrixFormatter.h"
#include "Parser/MatrixParser.h"


#include "/usr/local/include/gmp.h"
#include "mpreal.h"

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Matrix final {
	public:
		Matrix(){
			element[0][0] = 1; element[0][1] = 0; element[0][2] = 0;
			element[1][0] = 0; element[1][1] = 1; element[1][2] = 0;
			element[2][0] = 0; element[2][1] = 0; element[2][2] = 1;
		}
		
		Matrix(const T & m00, const T & m01,const T & m02,
					 const T & m10, const T & m11,const T & m12,
					 const T & m20, const T & m21,const T & m22){
			element[0][0] = m00; element[0][1] = m01; element[0][2] = m02;
			element[1][0] = m10; element[1][1] = m11; element[1][2] = m12;
			element[2][0] = m20; element[2][1] = m21; element[2][2] = m22;
		}
		
		~Matrix() {}
		
		static Matrix IdentityMatrix() { static Matrix matrix; return matrix; }
		
		void Element(const int & i, const int & j , const T & c){ element[i][j] = c; }
		
		T Element(const int & i, const int & j) const { return element[i][j]; }
		
		Vector<T> Line(const int & i) const{
			return Vector<T>(element[i][0],element[i][1],element[i][2]);
		}
		
		Vector<T> Column(const int & i) const {
			return Vector<T>(element[0][i],element[1][i],element[2][i]);
		}
		
		void Line(const int & i, const Vector<T> & l){
			element[i][0] = l.ComponantX();
			element[i][1] = l.ComponantY();
			element[i][2] = l.ComponantZ();
		}
		
		void Column(const int & i, const Vector<T> & c){
			element[0][i] = c.ComponantX();
			element[1][i] = c.ComponantY();
			element[2][i] = c.ComponantZ();
		}
		
		
		T Determinant() const{
			return element[0][0]*element[1][1]*element[2][2] + element[0][1]*element[1][2]*element[2][0] + element[0][2]*element[1][0]*element[2][1]
			- element[0][2]*element[1][1]*element[2][0] - element[1][2]*element[2][1]*element[0][0] - element[2][2]*element[1][0]*element[0][1];
		}
		
		Matrix<T> Product(const T & b) const{
			Matrix<T> a = *this;
			for(int i = 0 ; i < 3 ; i++)
				for(int j = 0 ; j < 3 ; j++)
					a.Element(i,j,a.Element(i,j)*b);
			return a;
		}
		
		Matrix<T> Div(const T & b) const{
			if(b == 0)
				throw(std::runtime_error("Dividing by 0 !"));
			Matrix<T> a = *this;
			for(int i = 0 ; i < 3 ; i++)
				for(int j = 0 ; j < 3 ; j++)
					a.Element(i,j,a.Element(i,j)/b);
			return a;
		}
		
		
		Matrix<T> MatrixTranspose() const{
			Matrix<T> a;
			Vector<T> raw;
			for(int i = 0 ; i < 3 ; i++){
				raw = Line(i);
				a.Column(i,raw);
			}
			return a;
		}
		
		Matrix<T> MatrixAdjoint() const{
			Matrix<T> a;
			a.Element(0,0,(element[1][1]*element[2][2]-element[1][2]*element[2][1]));
			a.Element(0,1,-(element[1][0]*element[2][2]-element[1][2]*element[2][0]));
			a.Element(0,2,(element[1][0]*element[2][1]-element[1][1]*element[2][0]));
			
			a.Element(1,0,-(element[0][1]*element[2][2]-element[0][2]*element[2][1]));
			a.Element(1,1,(element[0][0]*element[2][2]-element[0][2]*element[2][0]));
			a.Element(1,2,-(element[0][0]*element[2][1]-element[0][1]*element[2][0]));
			
			a.Element(2,0,(element[0][1]*element[1][2]-element[0][2]*element[1][1]));
			a.Element(2,1,-(element[0][0]*element[1][2]-element[0][2]*element[1][0]));
			a.Element(2,2,(element[0][0]*element[1][1]-element[1][0]*element[0][1]));
			return a;
		}
		
		Matrix<T> MatrixInverse() const{
			Matrix<T> a;
			T det = Determinant();
			if(det != 0){
				a = MatrixAdjoint().MatrixTranspose();
			}
			a /= det;
			return a;
		}
		
		
		Matrix<T> operator*(const T &b) const{
			Matrix<T> a = *this;
			a.Product(b);
			return a;
		}
		
		Vector<T> operator*(const Vector<T> &b) const{
			Vector<T> a;
			a.ComponantX(element[0][0]*b.ComponantX() + element[0][1]*b.ComponantY() + element[0][2]*b.ComponantZ());
			a.ComponantY(element[1][0]*b.ComponantX() + element[1][1]*b.ComponantY() + element[1][2]*b.ComponantZ());
			a.ComponantZ(element[2][0]*b.ComponantX() + element[2][1]*b.ComponantY() + element[2][2]*b.ComponantZ());
			return a;
		}
		
		Matrix<T> operator/(const T &b) const{
			Matrix<T> a = *this;
			a.Div(b);
			return a;
		}
		
		
		Matrix<T>& operator*=(const T & a){
			for(int i = 0 ; i < 3 ; i++)
				for(int j = 0 ; j < 3 ; j++)
					element[i][j]*=a;
			return *this;
		}
		
		Matrix<T>& operator/=(const T & a){
			if(a == 0)
				throw(std::runtime_error("Dividing by 0 !"));
			for(int i = 0 ; i < 3 ; i++)
				for(int j = 0 ; j < 3 ; j++)
					element[i][j]/=a;
			return *this;
		}
		
		
		std::string Format(MatrixFormatter<T>* matrixFormatter) {matrixFormatter->Format(*this);}
		
		void Parse(MatrixParser<T>* matrixParser, const std::string& str) {*this = matrixParser->Parse(str);}
				
	private:
		T element[3][3];
	};
	
}


template<class T> inline std::ostream & operator << (std::ostream & out, const GeometricalSpaceObjects::Matrix<T>& a){
	out << std::scientific << std::setprecision(16);
	for(int i = 0 ; i < 3 ; ++i){
		for(int j = 0 ; j < 3 ; ++j){
			out << a.Element(i,j);
			if(j != 2)
				out << "\t";
		}
		if(i != 2)
			out << "\n";
	}
	return out;
}

template<class T> inline std::istream & operator >> (std::istream & in, GeometricalSpaceObjects::Matrix<T> & a){
	T c;
	for(int i = 0 ; i < 3 ; ++i){
		for(int j = 0 ; j < 3 ; ++j){
			in >> c;
			a.Element(i,j,c);
		}
	}
	return in;
}

