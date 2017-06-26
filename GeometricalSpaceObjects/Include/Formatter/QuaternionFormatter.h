#pragma once

#include <string>
#include <sstream>

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Quaternion;
	
	template<class T>
	class QuaternionFormatter{
	public:
		virtual ~QuaternionFormatter() {}
		virtual std::string Format(const Quaternion<T>& quaternion) = 0;
		
	};
	
	template<class T>
	class LuGaQuaternionFormatter: public QuaternionFormatter<T>{
	public:
		LuGaQuaternionFormatter(){}
		~LuGaQuaternionFormatter(){}
		
		virtual std::string Format(const Quaternion<T>& quaternion) {
			std::stringstream sstr;
		  sstr.precision(15);
			sstr << std::scientific  << quaternion.ComponantReal() << "\t" << quaternion.ComponantI() << "\t" << quaternion.ComponantJ() << "\t" << quaternion.ComponantK();
			return sstr.str();
		}
	};

}
