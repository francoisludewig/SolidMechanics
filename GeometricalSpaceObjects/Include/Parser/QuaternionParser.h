#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Quaternion;
	
	template<class T>
	class QuaternionParser{
	public:
		virtual ~QuaternionParser() {}
		virtual Quaternion<T> Parse(const std::string& ) = 0;
		
	};
	
	template<class T>
	class LuGaQuaternionParser: public QuaternionParser<T>{
	public:
		LuGaQuaternionParser(){}
		~LuGaQuaternionParser(){}
		
		virtual Quaternion<T> Parse(const std::string& str) {
			std::stringstream sstr(str);
			T real = 0, i = 0, j = 0, k = 0;
			sstr >> real >> i >> i >> k;
			return Quaternion<T>(real, i, j, k);
		}
	};

}
