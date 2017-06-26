#pragma once 

#include <string>

#include "PointFormatter.h"
#include "QuaternionFormatter.h"

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Basis;
	
	template<class T>
	class BasisFormatter{
	public:
		virtual ~BasisFormatter() {}
		virtual std::string Format(const Basis<T>& basis) = 0;
		
	};

	template<class T>
	class LuGaBasisFormatter: public BasisFormatter<T>{
	public:
		LuGaBasisFormatter(){}
		~LuGaBasisFormatter(){}
		
		virtual std::string Format(const Basis<T>& basis) {
			auto origin = LuGaPointFormatter<T>().Format(basis.Origin());
			auto orientation = LuGaQuaternionFormatter<T>().Format(basis.Orientation());
			return origin + "\n" + orientation;
		}
	};
	
}
