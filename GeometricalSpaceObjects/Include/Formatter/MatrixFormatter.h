#pragma once

namespace GeometricalSpaceObjects {
	
	template<class T>
	class Matrix;
	
	template<class T>
	class MatrixFormatter{
	public:
		virtual ~MatrixFormatter() {}
		virtual std::string Format(const Matrix<T>& matrix) = 0;
		
	};
	
	template<class T>
	class LuGaMatrixFormatter: public MatrixFormatter<T>{
	public:
		LuGaMatrixFormatter(){}
		~LuGaMatrixFormatter(){}
		
		virtual std::string Format(const Matrix<T>& matrix) {;
			std::stringstream sstr;
			sstr.precision(15);
			sstr << std::scientific  << matrix.Element(0, 0) << "\t" << matrix.Element(0, 1) << "\t" << matrix.Element(0, 2) << "\n"
			<< matrix.Element(1, 0) << "\t" << matrix.Element(1, 1) << "\t" << matrix.Element(1, 2) << "\n"
			<< matrix.Element(2, 0) << "\t" << matrix.Element(2, 1) << "\t" << matrix.Element(2, 2);
			return sstr.str();
		}
	};

}
