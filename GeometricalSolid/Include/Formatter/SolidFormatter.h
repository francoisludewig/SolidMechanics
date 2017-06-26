#pragma once

#include "Solid.h"
#include <Formatter/BasisFormatter.h>
#include <Formatter/VectorFormatter.H>

using GeometricalSpaceObjects::LuGaBasisFormatter;
using GeometricalSpaceObjects::LuGaVectorFormatter;

namespace GeometricalSolid {
	
	class SolidFormatter{
	public:
		virtual ~SolidFormatter() {}
		virtual std::string Format(const Solid& solid) = 0;
		
	};
	
	class LuGaSolidFormatter: public SolidFormatter{
	public:
		LuGaSolidFormatter(){}
		~LuGaSolidFormatter(){}
		
		virtual std::string Format(const Solid& solid) {
			std::stringstream sstr;
			sstr.precision(15);
			sstr << basisFormatter.Format(solid.Basis()) << "\n";
			sstr << vectorFormatter.Format(solid.Force()) << "\n";
			sstr << vectorFormatter.Format(solid.Momentum());
			return sstr.str();
		}
	private:
		LuGaBasisFormatter<double> basisFormatter;
		LuGaVectorFormatter<double> vectorFormatter;
	};
	
}
