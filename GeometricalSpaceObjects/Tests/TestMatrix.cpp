#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

#include "Matrix.h"
#include "Precision.h"

using namespace std;
using namespace GeometricalSpaceObjects;

#define Matrix Matrix<Type>
#define Vector Vector<Type>


class MatrixTest : public ::testing::Test {
protected:
public:
	Matrix matrix;
	
	virtual void SetUp() {		
#ifndef DOUBLE_PRECISON
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(50));
#endif
		
		matrix.Element(0,0,pi/11);
		matrix.Element(0,1,pi);
		matrix.Element(0,2,pi/17);
		
		matrix.Element(1,0,pi/3);
		matrix.Element(1,1,2*pi);
		matrix.Element(1,2,pi);
		
		matrix.Element(2,0,pi);
		matrix.Element(2,1,pi/5);
		matrix.Element(2,2,pi/7);
	}
	
	virtual void TearDown() {}
};

TEST_F(MatrixTest,Constructor){
	Matrix a;
	
	EXPECT_MPREAL_EQ(1,a.Element(0,0));
	EXPECT_MPREAL_EQ(0,a.Element(0,1));
	EXPECT_MPREAL_EQ(0,a.Element(0,2));
	
	EXPECT_MPREAL_EQ(0,a.Element(1,0));
	EXPECT_MPREAL_EQ(1,a.Element(1,1));
	EXPECT_MPREAL_EQ(0,a.Element(1,2));
	
	EXPECT_MPREAL_EQ(0,a.Element(2,0));
	EXPECT_MPREAL_EQ(0,a.Element(2,1));
	EXPECT_MPREAL_EQ(1,a.Element(2,2));
}

TEST_F(MatrixTest,det){
	Matrix a;
	
	a.Element(0,0,-3);
	a.Element(0,1,5);
	a.Element(0,2,6);
	
	a.Element(1,0,-1);
	a.Element(1,1,2);
	a.Element(1,2,2);
	
	a.Element(2,0,1);
	a.Element(2,1,-1);
	a.Element(2,2,-1);
	EXPECT_MPREAL_EQ(-1,a.Determinant());
}

TEST_F(MatrixTest,Transpose){
	Matrix a,b;
	
	a.Element(0,0,-3);
	a.Element(0,1,5);
	a.Element(0,2,6);
	
	a.Element(1,0,-1);
	a.Element(1,1,2);
	a.Element(1,2,2);
	
	a.Element(2,0,1);
	a.Element(2,1,-1);
	a.Element(2,2,-1);
	
	b = a.MatrixTranspose();
	
	EXPECT_MPREAL_EQ(-3,b.Element(0,0));
	EXPECT_MPREAL_EQ(-1,b.Element(0,1));
	EXPECT_MPREAL_EQ(1,b.Element(0,2));
	
	EXPECT_MPREAL_EQ(5,b.Element(1,0));
	EXPECT_MPREAL_EQ(2,b.Element(1,1));
	EXPECT_MPREAL_EQ(-1,b.Element(1,2));
	
	EXPECT_MPREAL_EQ(6,b.Element(2,0));
	EXPECT_MPREAL_EQ(2,b.Element(2,1));
	EXPECT_MPREAL_EQ(-1,b.Element(2,2));
}


TEST_F(MatrixTest,Inverse){
	Matrix a,b;
	
	a.Element(0,0,-3);
	a.Element(0,1,5);
	a.Element(0,2,6);
	
	a.Element(1,0,-1);
	a.Element(1,1,2);
	a.Element(1,2,2);
	
	a.Element(2,0,1);
	a.Element(2,1,-1);
	a.Element(2,2,-1);
	
	b = a.MatrixInverse();
	
	EXPECT_MPREAL_EQ(0,b.Element(0,0));
	EXPECT_MPREAL_EQ(1,b.Element(0,1));
	EXPECT_MPREAL_EQ(2,b.Element(0,2));
	
	EXPECT_MPREAL_EQ(-1,b.Element(1,0));
	EXPECT_MPREAL_EQ(3,b.Element(1,1));
	EXPECT_MPREAL_EQ(0,b.Element(1,2));
	
	EXPECT_MPREAL_EQ(1,b.Element(2,0));
	EXPECT_MPREAL_EQ(-2,b.Element(2,1));
	EXPECT_MPREAL_EQ(1,b.Element(2,2));
}


TEST_F(MatrixTest,VectorProductIdentity){
	Matrix a;
	Vector v(pi,pi/2,pi/3);
	Vector b = a*v;
	
	EXPECT_MPREAL_EQ(v.ComponantX(),b.ComponantX());
	EXPECT_MPREAL_EQ(v.ComponantY(),b.ComponantY());
	EXPECT_MPREAL_EQ(v.ComponantZ(),b.ComponantZ());
}

TEST_F(MatrixTest,VectorProduct){
	Matrix a(1,2,3,5,4,6,7,-4,-3);
	Vector v(pi,pi/2,pi/3);
	Vector b = a*v;
	
	EXPECT_MPREAL_EQ(3*pi,b.ComponantX());
	EXPECT_MPREAL_EQ(9*pi,b.ComponantY());
	EXPECT_MPREAL_EQ(4*pi,b.ComponantZ());
}

TEST_F(MatrixTest,Formatter){
	const std::string expected("2.8559933214452665e-01\t3.1415926535897931e+00\t1.8479956785822313e-01\n1.0471975511965976e+00\t6.2831853071795862e+00\t3.1415926535897931e+00\n3.1415926535897931e+00\t6.2831853071795862e-01\t4.4879895051282759e-01");
	std::stringstream sstr;
	sstr << matrix;
	auto text = sstr.str();
	std::cout << text << std::endl;

  EXPECT_TRUE(expected == text);
}

TEST_F(MatrixTest,Parser){
		const std::string expected("2.8559933214452665e-01\t3.1415926535897931e+00\t1.8479956785822313e-01\n1.0471975511965976e+00\t6.2831853071795862e+00\t3.1415926535897931e+00\n3.1415926535897931e+00\t6.2831853071795862e-01\t4.4879895051282759e-01");
	std::stringstream sstr(expected);

	Matrix mat;
	sstr >> mat;
	
	EXPECT_MPREAL_EQ(pi/11, mat.Element(0,0));
	EXPECT_MPREAL_EQ(pi, mat.Element(0,1));
	EXPECT_MPREAL_EQ(pi/17, mat.Element(0,2));
	
	EXPECT_MPREAL_EQ(pi/3, mat.Element(1,0));
	EXPECT_MPREAL_EQ(2*pi, mat.Element(1,1));
	EXPECT_MPREAL_EQ(pi, mat.Element(1,2));
	
	EXPECT_MPREAL_EQ(pi, mat.Element(2,0));
	EXPECT_MPREAL_EQ(pi/5, mat.Element(2,1));
	EXPECT_MPREAL_EQ(pi/7, mat.Element(2,2));
}

TEST_F(MatrixTest,FormatterParser){
	const std::string expected("2.855993321445267e-01\t3.141592653589793e+00\t1.847995678582231e-01\n1.047197551196598e+00\t6.283185307179586e+00\t3.141592653589793e+00\n3.141592653589793e+00\t6.283185307179586e-01\t4.487989505128276e-01");
	std::stringstream sstr;
	sstr << matrix;
	auto text = sstr.str();
	std::stringstream sstr2(text);
	
	Matrix mat;
	sstr2 >> mat;
	
	EXPECT_MPREAL_EQ(matrix.Element(0,0), mat.Element(0,0));
	EXPECT_MPREAL_EQ(matrix.Element(0,1), mat.Element(0,1));
	EXPECT_MPREAL_EQ(matrix.Element(0,2), mat.Element(0,2));
	
	EXPECT_MPREAL_EQ(matrix.Element(1,0), mat.Element(1,0));
	EXPECT_MPREAL_EQ(matrix.Element(1,1), mat.Element(1,1));
	EXPECT_MPREAL_EQ(matrix.Element(1,2), mat.Element(1,2));
	
	EXPECT_MPREAL_EQ(matrix.Element(2,0), mat.Element(2,0));
	EXPECT_MPREAL_EQ(matrix.Element(2,1), mat.Element(2,1));
	EXPECT_MPREAL_EQ(matrix.Element(2,2), mat.Element(2,2));
}

TEST_F(MatrixTest,IO_Operator){
	Matrix a;
	Matrix b;
	
	for(int i = 0 ; i < 3 ; i++){
		for(int j = 0 ; j < 3 ; j++){
			a.Element(i,j,pi*(i+1)/(j+1));
		}
	}
	ofstream fichierOut("testMatrix3x3.txt", ios::out | ios::trunc);
	if(fichierOut){
		fichierOut << a;
		fichierOut.close();
	}
	
	ifstream fichierIn("testMatrix3x3.txt", ios::in);
	if(fichierIn){
		fichierIn >> b;
		fichierIn.close();
	}
	
	for(int i = 0 ; i < 3 ; i++){
		for(int j = 0 ; j < 3 ; j++){
			EXPECT_MPREAL_EQ(a.Element(i,j),b.Element(i,j));
		}
	}
	remove("testMatrix3x3.txt");
}
