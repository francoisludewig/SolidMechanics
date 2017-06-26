#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Point;
	
	template<class T>
	class PointFormatter{
	public:
		virtual ~PointFormatter() {}
		virtual std::string Format(const Point<T>& point) = 0;
	
	};
	
	template<class T>
	class LuGaPointFormatter: public PointFormatter<T>{
	public:
		LuGaPointFormatter(){}
		~LuGaPointFormatter(){}
		
		virtual std::string Format(const Point<T>& point) {
			std::stringstream sstr;
			sstr.precision(15);
			sstr << std::scientific  << point.CoordinateX() << "\t" << point.CoordinateY() << "\t" << point.CoordinateZ();
			return sstr.str();
		}
	};
	
}
