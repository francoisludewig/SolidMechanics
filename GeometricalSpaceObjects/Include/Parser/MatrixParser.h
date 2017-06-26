#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Matrix;
	
	template<class T>
	class MatrixParser{
	public:
		virtual ~MatrixParser() {}
		virtual Matrix<T> Parse(const std::string& ) = 0;
		
	};
	
	template<class T>
	class LuGaMatrixParser: public MatrixParser<T>{
	public:
		LuGaMatrixParser(){}
		~LuGaMatrixParser(){}
		
		virtual Matrix<T> Parse(const std::string& str) {
			std::stringstream sstr(str);
			T element;
			Matrix<T> matrix;
			for(int i = 0 ; i < 3 ; ++i) {
				for(int j = 0 ; j < 3 ; ++j) {
					sstr >> element;
					matrix.Element(i, j, element);
				}
			}
			return std::move(matrix);
		}
	};

}
