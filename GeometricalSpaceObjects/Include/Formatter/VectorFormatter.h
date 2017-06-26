#pragma once

#include <sstream>

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Vector;
	
	template<class T>
	class VectorFormatter{
	public:
		virtual ~VectorFormatter() {}
		virtual std::string Format(const Vector<T>& vector) = 0;
		
	};
	
	template<class T>
	class LuGaVectorFormatter: public VectorFormatter<T>{
	public:
		LuGaVectorFormatter(){}
		~LuGaVectorFormatter(){}
		
		virtual std::string Format(const Vector<T>& vector) {
			std::stringstream sstr;
			sstr.precision(15);
			sstr << std::scientific  << vector.ComponantX() << "\t" << vector.ComponantY() << "\t" << vector.ComponantZ();
			return sstr.str();
		}
	};

}
