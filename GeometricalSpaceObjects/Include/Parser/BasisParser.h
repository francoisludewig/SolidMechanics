#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Basis;
	
	template<class T>
	class BasisParser{
	public:
		virtual ~BasisParser() {}
		virtual Basis<T> Parse(const std::string&) = 0;
		
	};
	
	template<class T>
	class LuGaBasisParser: public BasisParser<T>{
	public:
		LuGaBasisParser(){}
		~LuGaBasisParser(){}
		
		virtual Basis<T> Parse(const std::string& str) {
			std::stringstream sstr(str);
			Basis<T> basis;
			T x = 0,y = 0,z = 0;
			sstr >> x >> y >> z;
			basis.Origin(Point<T>(x, y, z));
			T real = 0, i = 0, j = 0, k = 0;
			sstr >> real >> i >> i >> k;
			basis.Orientation(Quaternion<T>(real, i, j, k));
			return std::move(basis);
		}
	};
	
}
