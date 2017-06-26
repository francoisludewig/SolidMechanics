#pragma once

#include <iostream>
#include <iomanip>

#include "Formatter/PointFormatter.h"
#include "Parser/PointParser.h"

#include "/usr/local/include/gmp.h"
#include "mpreal.h"

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Vector;
	
	template<class T>
	class Point final {
	public:
		Point():coordinateX(0),coordinateY(0),coordinateZ(0) {}

		Point(const T &x, const T &y, const T &z):coordinateX(x),coordinateY(y),coordinateZ(z) {}

		~Point() {}
		T CoordinateX() const{return this->coordinateX;}
		T CoordinateY() const{return this->coordinateY;}
		T CoordinateZ() const{return this->coordinateZ;}
		void SetCoordinates(const T & x, const T & y, const T & z){
			this->coordinateX = x;
			this->coordinateY = y;
			this->coordinateZ = z;
		}

		void Translate(const Vector<T> & a){
			coordinateX += a.ComponantX();
			coordinateY += a.ComponantY();
			coordinateZ += a.ComponantZ();
		}

		
		Point operator+(const Vector<T> &b) const{
			return std::move(Point<T> (b.ComponantX()+coordinateX,b.ComponantY()+coordinateY,b.ComponantZ()+coordinateZ));
		}

		Vector<T> operator-(const Point &b) const{
			return std::move(Vector<T>(coordinateX-b.CoordinateX(),coordinateY-b.CoordinateY(),coordinateZ-b.CoordinateZ()));
		}

		void operator+=(const Vector<T>& a) { Translate(a); }
		
		void operator-=(const Vector<T>& a) { Translate(-a); }
		
		
		std::string Format(PointFormatter<T>* pointFormatter) {pointFormatter->Format(*this);}
		
		void Parse(PointParser<T>* pointParser, const std::string& str) {*this = pointParser->Parse(str);}
		
		
	private:
		T coordinateX,coordinateY,coordinateZ;
	};
	
		template<class T>
	bool operator== (const Point<T>& pt1, const Point<T>& pt2){
		if(fabs(pt1.CoordinateX()-pt2.CoordinateX()) > std::numeric_limits<T>::epsilon())
			return false;
		
		if(fabs(pt1.CoordinateY()-pt2.CoordinateY()) > std::numeric_limits<T>::epsilon())
			return false;
		
		if(fabs(pt1.CoordinateZ()-pt2.CoordinateZ()) > std::numeric_limits<T>::epsilon())
			return false;
		
		return true;
	}

	template<class T>
	bool operator!= (const Point<T>& pt1, const Point<T>& pt2){
		return !(pt1==pt2);
	}
	
}

template<class T> inline std::ostream & operator << (std::ostream & out, const GeometricalSpaceObjects::Point<T> & a){
	out << std::scientific << std::setprecision(16);
	out << a.CoordinateX() << "\t" << a.CoordinateY() << "\t" << a.CoordinateZ();
	return out;
}

template<class T> inline std::istream & operator >> (std::istream & in, GeometricalSpaceObjects::Point<T> & a){
	T x,y,z;
	in >> x >> y >> z;
	a.SetCoordinates(x,y,z);
	return in;
}

