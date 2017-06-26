#pragma once 

#include "Solid.h"
#include <string>
#include <sstream>
#include <Parser/BasisParser.h>
#include <Parser/VectorParser.h>


namespace GeometricalSolid {
	
	class Solid;
	
	class SolidParser{
	public:
		virtual ~SolidParser() {}
		virtual Solid Parse(const std::string&) = 0;
		
	};
	/*
	class LuGaSolidParser: public SolidParser{
	public:
		LuGaSolidParser(){}
		~LuGaSolidParser(){}
		
		virtual Solid Parse(const std::string&) {
			Solid s;		
			return s;
		}
	};
	 */
	
}
