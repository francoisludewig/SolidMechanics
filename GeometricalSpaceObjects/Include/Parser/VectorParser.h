#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Vector;
	
	template<class T>
	class VectorParser{
	public:
		virtual ~VectorParser() {}
		virtual Vector<T> Parse(const std::string&) = 0;
		
	};
	
	template<class T>
	class LuGaVectorParser: public VectorParser<T>{
	public:
		LuGaVectorParser(){}
		~LuGaVectorParser(){}
		
		virtual Vector<T> Parse(const std::string& str) {
			std::stringstream sstr(str);
			T x = 0,y = 0,z = 0;
			sstr >> x >> y >> z;
			return Vector<T>(x, y, z);
		}
	};

	
	
}
