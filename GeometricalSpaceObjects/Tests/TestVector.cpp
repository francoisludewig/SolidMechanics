#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

#include "Vector.h"
#include "Precision.h"

using namespace std;
using namespace GeometricalSpaceObjects;

#define Vector Vector<Type>

class VectorTest : public ::testing::Test {
protected:
	VectorTest() {
		srand(time(NULL));
	}
	
	virtual void SetUp() {
#ifndef DOUBLE_PRECISON
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(50));
#endif
	}
	
	virtual void TearDown() {}
};


TEST_F(VectorTest,Constructor){
	Vector a;
	EXPECT_MPREAL_EQ(0,a.ComponantX());
	EXPECT_MPREAL_EQ(0,a.ComponantY());
	EXPECT_MPREAL_EQ(0,a.ComponantZ());
}

TEST_F(VectorTest,ConstructorValue){
	Vector a(1,2,3);
	EXPECT_MPREAL_EQ(1,a.ComponantX());
	EXPECT_MPREAL_EQ(2,a.ComponantY());
	EXPECT_MPREAL_EQ(3,a.ComponantZ());
	
	Vector b(pi,2*pi,pi/2);
	EXPECT_MPREAL_EQ(pi,b.ComponantX());
	EXPECT_MPREAL_EQ(2*pi,b.ComponantY());
	EXPECT_MPREAL_EQ(pi/2,b.ComponantZ());
}

TEST_F(VectorTest,Norme){
	Vector a(1,2,3);
	EXPECT_MPREAL_EQ(sqrt(14.),a.Norme());
}


TEST_F(VectorTest,Normalize){
	Vector a(1,2,3),b = a;
	a.Normalize();
	EXPECT_MPREAL_EQ(1,a.Norme());
}

TEST_F(VectorTest,ScalarProduct_WhenParallel){
	Vector a(1,2,3),b = a;
	a.Normalize();
	Vector pv = a^b;
	EXPECT_TRUE(pv == Vector(0,0,0));
}


TEST_F(VectorTest,ScalarProduct){
	Vector a(1,0,0);
	Vector b(0,1,0);
	Vector c(0,0,1);
	
	EXPECT_MPREAL_EQ(0,a*b);
	EXPECT_MPREAL_EQ(0,a*c);
	EXPECT_MPREAL_EQ(0,b*c);
	EXPECT_MPREAL_EQ(1,a*a);
	EXPECT_MPREAL_EQ(1,b*b);
	EXPECT_MPREAL_EQ(1,c*c);
}


TEST_F(VectorTest,VectorialProduct){
	Vector a(1,0,0);
	Vector b(0,1,0);
	Vector c(0,0,1);
	
	Vector expectedC = a^b;
	Vector expectedA = b^c;
	Vector expectedB = c^a;
	
	EXPECT_TRUE(c == expectedC);
	EXPECT_TRUE(a == expectedA);
	EXPECT_TRUE(b == expectedB);
}


TEST_F(VectorTest,Product){
	Vector a(1,-2,3);
	Vector b(6,9,-7);
	
	EXPECT_MPREAL_EQ(0,a*(a^b));
	EXPECT_MPREAL_EQ(0,b*(a^b));
	EXPECT_MPREAL_EQ(0,a*(b^a));
	EXPECT_MPREAL_EQ(0,b*(b^a));
	
	a.SetComponants(0,0,0);
	b.SetComponants(1,3,1);
	
	EXPECT_MPREAL_EQ(0,a*(a^b));
	EXPECT_MPREAL_EQ(0,b*(a^b));
	EXPECT_MPREAL_EQ(0,a*(b^a));
	EXPECT_MPREAL_EQ(0,b*(b^a));
}


TEST_F(VectorTest,ProductUnit){
	Vector a(pi,2*pi,pi/2);
	Vector b(a.ComponantY()*a.ComponantZ(),a.ComponantX()*a.ComponantZ(),-2*a.ComponantX()*a.ComponantY());
	a.Normalize();
	b.Normalize();
	Vector c = a^b;
	
	EXPECT_MPREAL_EQ(0,b*a);
	EXPECT_MPREAL_EQ(0,c*a);
	EXPECT_MPREAL_EQ(0,c*b);
	EXPECT_MPREAL_EQ(1,c.Norme());
}


TEST_F(VectorTest,MultiplicationEqual){
	Vector a(5,17,-57);
	double b = -13;
	Vector c = a*b;
	a *=b;
	
	EXPECT_MPREAL_EQ(a.ComponantX(),c.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY(),c.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ(),c.ComponantZ());
}


TEST_F(VectorTest,IO_Operator){
	Vector a(pi,2*pi,pi/2);
	Vector b;
	
	ofstream fichierOut("testVector3D.txt", ios::out | ios::trunc);
	if(fichierOut){
		fichierOut << a;
		fichierOut.close();
	}
	
	ifstream fichierIn("testVector3D.txt", ios::in);
	if(fichierIn){
		fichierIn >> b;
		fichierIn.close();
	}
	
	EXPECT_MPREAL_EQ(a.ComponantX(),b.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY(),b.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ(),b.ComponantZ());
	
	remove("testVector3D.txt");
}

TEST_F(VectorTest,Sum){
	Vector a(pi,2*pi,pi/2);
	Vector b(7*pi,-2*pi,-4*pi/2);
	Vector c = a+b;
	
	EXPECT_MPREAL_EQ(a.ComponantX()+b.ComponantX(),c.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY()+b.ComponantY(),c.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ()+b.ComponantZ(),c.ComponantZ());
}


TEST_F(VectorTest,Diff){
	Vector a(pi,2*pi,pi/2);
	Vector b(7*pi,2*pi,-4*pi/2);
	Vector c = a-b;
	
	EXPECT_MPREAL_EQ(a.ComponantX()-b.ComponantX(),c.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY()-b.ComponantY(),c.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ()-b.ComponantZ(),c.ComponantZ());
}

TEST_F(VectorTest,Div){
	Vector a(pi,2*pi,pi/2);
	double b = sqrt(2.);
	Vector c = a/b;
	
	EXPECT_MPREAL_EQ(a.ComponantX()/b,c.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY()/b,c.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ()/b,c.ComponantZ());
}


TEST_F(VectorTest,DivEqual){
	Vector a(pi,2*pi,pi/2);
	double b = sqrt(2.);
	Vector c = a/b;
	a /=b;
	EXPECT_MPREAL_EQ(a.ComponantX(),c.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY(),c.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ(),c.ComponantZ());
}


TEST_F(VectorTest,LinearCombination){
	Vector a(1,0,0);
	Vector b(0,1,0);
	Vector c(0,0,1);
	
	Vector d = static_cast<Type>(sqrt(3.))*a+static_cast<Type>(sqrt(5.))*b-static_cast<Type>(sqrt(7.))*c;
	EXPECT_MPREAL_EQ(sqrt(3.0),d.ComponantX());
	EXPECT_MPREAL_EQ(sqrt(5.0),d.ComponantY());
	EXPECT_MPREAL_EQ(-sqrt(7.0),d.ComponantZ());
}


TEST_F(VectorTest,increment){
	Vector a(pi,0,0);
	a += static_cast<Type>(2)*a;
	EXPECT_MPREAL_EQ(3*pi,a.ComponantX());
}

TEST_F(VectorTest,decrement){
	Vector a(pi,0,0);
	a -= static_cast<Type>(2)*a;
	EXPECT_MPREAL_EQ(-pi,a.ComponantX());
}

TEST_F(VectorTest,ScalarDecomposition){
	Vector a(1,0,0);
	Vector b(0,1,0);
	Vector c(0,0,1);
	Vector d(5*pi,-7*pi,13*pi);
	
	EXPECT_MPREAL_EQ(5*pi,d*a);
	EXPECT_MPREAL_EQ(-7*pi,d*b);
	EXPECT_MPREAL_EQ(13*pi,d*c);
}

TEST_F(VectorTest,VectorialRandom){
	Vector a(pi,2*pi,pi/2);
	Vector b(pi/6,pi/3,pi/4);
	
	Vector pvExpected = a^b;
	
	EXPECT_MPREAL_EQ(a.ComponantY()*b.ComponantZ()-a.ComponantZ()*b.ComponantY(),pvExpected.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantZ()*b.ComponantX()-a.ComponantX()*b.ComponantZ(),pvExpected.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantX()*b.ComponantY()-a.ComponantY()*b.ComponantX(),pvExpected.ComponantZ());
}

TEST_F(VectorTest,Formatter){
	Vector a(pi,2*pi,pi/2);
	std::stringstream sstr;
	sstr << a;
	auto text = sstr.str();
	EXPECT_TRUE(std::string("3.1415926535897931e+00	6.2831853071795862e+00	1.5707963267948966e+00") == text);
}

TEST_F(VectorTest,Parser){
	std::stringstream sstr(std::string("3.1415926535897931e+00	6.2831853071795862e+00	1.5707963267948966e+00"));
	Vector vect;
	sstr >> vect;
	EXPECT_MPREAL_EQ(pi, vect.ComponantX());
	EXPECT_MPREAL_EQ(2*pi, vect.ComponantY());
	EXPECT_MPREAL_EQ(pi/2, vect.ComponantZ());
}

TEST_F(VectorTest,FormatterParser){
	Vector a(pi,2*pi,pi/2);
	std::stringstream sstr;
	sstr << a;
	auto text = sstr.str();
	
	std::stringstream sstr2(text);
	Vector vect;
	sstr2 >> vect;

	EXPECT_MPREAL_EQ(a.ComponantX(), vect.ComponantX());
	EXPECT_MPREAL_EQ(a.ComponantY(), vect.ComponantY());
	EXPECT_MPREAL_EQ(a.ComponantZ(), vect.ComponantZ());
}
 
