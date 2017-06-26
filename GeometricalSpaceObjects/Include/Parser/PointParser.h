#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Point;
	
	template<class T>
	class PointParser{
	public:
		virtual ~PointParser() {}
		virtual Point<T> Parse(const std::string& ) = 0;
		
	};
	
	template<class T>
	class LuGaPointParser: public PointParser<T>{
	public:
		LuGaPointParser(){}
		~LuGaPointParser(){}
		
		virtual Point<T> Parse(const std::string& str) {
			std::stringstream sstr(str);
			T x = 0,y = 0,z = 0;
			sstr >> x >> y >> z;
			return Point<T>(x, y, z);
		}
	};

}
