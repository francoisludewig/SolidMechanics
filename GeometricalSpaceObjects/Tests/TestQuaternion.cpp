#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

#include "Quaternion.h"
#include "Vector.h"
#include "Precision.h"

using namespace std;
using namespace GeometricalSpaceObjects;

#define Quaternion Quaternion<Type>

class QuaternionTest : public ::testing::Test {
protected:
	virtual void SetUp() {
#ifndef DOUBLE_PRECISON
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(50));
#endif
	}
	
	virtual void TearDown() {}
};


TEST_F(QuaternionTest,Constructor){
	Quaternion a;
	EXPECT_MPREAL_EQ(1,a.ComponantReal());
	EXPECT_MPREAL_EQ(0,a.ComponantI());
	EXPECT_MPREAL_EQ(0,a.ComponantJ());
	EXPECT_MPREAL_EQ(0,a.ComponantK());
}

TEST_F(QuaternionTest,ConstructorValue){
	Quaternion a(pi,2*pi,pi/2,pi/3);
	EXPECT_MPREAL_EQ(pi,a.ComponantReal());
	EXPECT_MPREAL_EQ(2*pi,a.ComponantI());
	EXPECT_MPREAL_EQ(pi/2,a.ComponantJ());
	EXPECT_MPREAL_EQ(pi/3,a.ComponantK());
}

TEST_F(QuaternionTest,ConstructorVector3D){
	Vector<Type> w(pi,0,0);
	Quaternion a(w);
	
	EXPECT_MPREAL_EQ(0,a.ComponantReal());
	EXPECT_MPREAL_EQ(1.,a.ComponantI());
	EXPECT_MPREAL_EQ(0.,a.ComponantJ());
	EXPECT_MPREAL_EQ(0.,a.ComponantK());
}


TEST_F(QuaternionTest,SetComponants){
	Quaternion a;
	a.SetComponants(0.5,0.5,0.5,0.5);
	EXPECT_MPREAL_EQ(0.5,a.ComponantReal());
	EXPECT_MPREAL_EQ(0.5,a.ComponantI());
	EXPECT_MPREAL_EQ(0.5,a.ComponantJ());
	EXPECT_MPREAL_EQ(0.5,a.ComponantK());
}


TEST_F(QuaternionTest,multiplication){
	Quaternion a,b,c;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	b.SetComponants(6*pi,pi/4.,pi/6,pi);
	c = a*b;
	
	EXPECT_MPREAL_EQ(c.ComponantReal(),a.ComponantReal()*b.ComponantReal() - a.ComponantI()*b.ComponantI() - a.ComponantJ()*b.ComponantJ() - a.ComponantK()*b.ComponantK());
	EXPECT_MPREAL_EQ(c.ComponantI(),a.ComponantReal()*b.ComponantI() + a.ComponantI()*b.ComponantReal() - a.ComponantJ()*b.ComponantK() + a.ComponantK()*b.ComponantJ());
	EXPECT_MPREAL_EQ(c.ComponantJ(),a.ComponantReal()*b.ComponantJ() + a.ComponantI()*b.ComponantK() + a.ComponantJ()*b.ComponantReal() - a.ComponantK()*b.ComponantI());
	EXPECT_MPREAL_EQ(c.ComponantK(),a.ComponantReal()*b.ComponantK() - a.ComponantI()*b.ComponantJ() + a.ComponantJ()*b.ComponantI() + a.ComponantK()*b.ComponantReal());
}


TEST_F(QuaternionTest,Sum){
	Quaternion a,b,c;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	b.SetComponants(6*pi,pi/4.,pi/6,pi);
	c = a+b;
	
	EXPECT_MPREAL_EQ(c.ComponantReal(),a.ComponantReal()+b.ComponantReal());
	EXPECT_MPREAL_EQ(c.ComponantI(),a.ComponantI()+b.ComponantI());
	EXPECT_MPREAL_EQ(c.ComponantJ(),a.ComponantJ()+b.ComponantJ());
	EXPECT_MPREAL_EQ(c.ComponantK(),a.ComponantK()+b.ComponantK());
}


TEST_F(QuaternionTest,Diff){
	Quaternion a,b,c;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	b.SetComponants(6*pi,pi/4.,pi/6,pi);
	c = a-b;
	
	EXPECT_MPREAL_EQ(c.ComponantReal(),a.ComponantReal()-b.ComponantReal());
	EXPECT_MPREAL_EQ(c.ComponantI(),a.ComponantI()-b.ComponantI());
	EXPECT_MPREAL_EQ(c.ComponantJ(),a.ComponantJ()-b.ComponantJ());
	EXPECT_MPREAL_EQ(c.ComponantK(),a.ComponantK()-b.ComponantK());
}



TEST_F(QuaternionTest,SumOperateurUnaire){
	Quaternion a,b;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	b.SetComponants(6*pi,pi/4.,pi/6,pi);
	a+=b;
	
	EXPECT_MPREAL_EQ(a.ComponantReal(),7*pi);
	EXPECT_MPREAL_EQ(a.ComponantI(),9*pi/4);
	EXPECT_MPREAL_EQ(a.ComponantJ(),2*pi/3);
	EXPECT_MPREAL_EQ(a.ComponantK(),4*pi/3);
}


TEST_F(QuaternionTest,DiffOperateurUnaire){
	Quaternion a,b;
	a.SetComponants(pi   , 2*pi , pi/2 , pi/3);
	b.SetComponants(6*pi , pi/4 , pi/6 , pi);
	a-=b;
	
	EXPECT_MPREAL_EQ(a.ComponantReal(),-5*pi);
	EXPECT_MPREAL_EQ(a.ComponantI(),7*pi/4);
	EXPECT_MPREAL_EQ(a.ComponantJ(),pi/3);
	EXPECT_MPREAL_EQ(a.ComponantK(),(-2/3)*pi);
}




TEST_F(QuaternionTest,SumConjugation){
	Quaternion a,c;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	c = a+~a;
	
	EXPECT_MPREAL_EQ(c.ComponantReal(),2*a.ComponantReal());
	EXPECT_MPREAL_EQ(c.ComponantI(),0);
	EXPECT_MPREAL_EQ(c.ComponantJ(),0);
	EXPECT_MPREAL_EQ(c.ComponantK(),0);
}


TEST_F(QuaternionTest,DiffConjugation){
	Quaternion a,c;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	c = a-~a;
	
	EXPECT_MPREAL_EQ(c.ComponantReal(),0);
	EXPECT_MPREAL_EQ(c.ComponantI(),2*a.ComponantI());
	EXPECT_MPREAL_EQ(c.ComponantJ(),2*a.ComponantJ());
	EXPECT_MPREAL_EQ(c.ComponantK(),2*a.ComponantK());
}



TEST_F(QuaternionTest,multiplicationOperator){
	Quaternion a;
	a.SetComponants(0.5,0.5,0.5,0.5);
	Quaternion b,c;
	
	c = a*b;
	EXPECT_MPREAL_EQ(c.ComponantReal(),a.ComponantReal());
	EXPECT_MPREAL_EQ(c.ComponantI(),a.ComponantI());
	EXPECT_MPREAL_EQ(c.ComponantJ(),a.ComponantJ());
	EXPECT_MPREAL_EQ(c.ComponantK(),a.ComponantK());
	
	c = b*a;
	EXPECT_MPREAL_EQ(c.ComponantReal(),a.ComponantReal());
	EXPECT_MPREAL_EQ(c.ComponantI(),a.ComponantI());
	EXPECT_MPREAL_EQ(c.ComponantJ(),a.ComponantJ());
	EXPECT_MPREAL_EQ(c.ComponantK(),a.ComponantK());
	
	a*=b;
	EXPECT_MPREAL_EQ(0.5,a.ComponantReal());
	EXPECT_MPREAL_EQ(0.5,a.ComponantI());
	EXPECT_MPREAL_EQ(0.5,a.ComponantJ());
	EXPECT_MPREAL_EQ(0.5,a.ComponantK());
}


TEST_F(QuaternionTest,conjuguation){
	Quaternion a;
	a.SetComponants(0.5,0.5,0.5,0.5);
	Quaternion b;
	b = ~a;
	
	EXPECT_MPREAL_EQ(b.ComponantReal(),a.ComponantReal());
	EXPECT_MPREAL_EQ(b.ComponantI(),-a.ComponantI());
	EXPECT_MPREAL_EQ(b.ComponantJ(),-a.ComponantJ());
	EXPECT_MPREAL_EQ(b.ComponantK(),-a.ComponantK());
}


TEST_F(QuaternionTest,conjuguationAndMultiplication){
	Quaternion a;
	a.SetComponants(0.5,0.5,0.5,0.5);
	Quaternion b,c;
	b = ~a;
	c = a*b;
	
	EXPECT_MPREAL_EQ(1,c.ComponantReal());
	EXPECT_MPREAL_EQ(0,c.ComponantI());
	EXPECT_MPREAL_EQ(0,c.ComponantJ());
	EXPECT_MPREAL_EQ(0,c.ComponantK());
}



TEST_F(QuaternionTest,IO_Operator){
	Quaternion a;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	Quaternion b;
	
	ofstream fichierOut("testUnitQuaternion.txt", ios::out | ios::trunc);
	if(fichierOut){
		fichierOut << a;
		fichierOut.close();
	}
	
	ifstream fichierIn("testUnitQuaternion.txt", ios::in);
	if(fichierIn){
		fichierIn >> b;
		fichierIn.close();
	}
	
	EXPECT_MPREAL_EQ(a.ComponantReal(),b.ComponantReal());
	EXPECT_MPREAL_EQ(a.ComponantI(),b.ComponantI());
	EXPECT_MPREAL_EQ(a.ComponantJ(),b.ComponantJ());
	EXPECT_MPREAL_EQ(a.ComponantK(),b.ComponantK());
	
	remove("testUnitQuaternion.txt");
}

TEST_F(QuaternionTest,Multiplications){
	Quaternion a;
	a.SetComponants(pi,2*pi,pi/2,pi/3);
	a.Normalize();
	Quaternion b;
	
	for(int i = 0 ; i < 1000000 ; ++i){
		b.SetComponants(static_cast<Type>((rand()%RAND_MAX)) / static_cast<Type>(RAND_MAX),
										static_cast<Type>((rand()%RAND_MAX)) / static_cast<Type>(RAND_MAX),
										static_cast<Type>((rand()%RAND_MAX)) / static_cast<Type>(RAND_MAX),
										static_cast<Type>((rand()%RAND_MAX)) / static_cast<Type>(RAND_MAX));
		b.Normalize();
		a *= b;
		EXPECT_MPREAL_EQ(1.,a.Norme());
	}
}


TEST_F(QuaternionTest,Format){
	Quaternion a(pi,2*pi,pi/2,pi/3);
	a.Normalize();
	std::stringstream sstr;
	sstr << a;
	auto text = sstr.str();
	EXPECT_TRUE(std::string("4.3188945044921667e-01	8.6377890089843334e-01	2.1594472522460834e-01	1.4396315014973887e-01") == text);
}

TEST_F(QuaternionTest,Parser){
	std::stringstream sstr(std::string("4.3188945044921667e-01	8.6377890089843334e-01	2.1594472522460834e-01	1.4396315014973887e-01"));
	Quaternion quad;
	sstr >> quad;
	
	Type norme = sqrt(pi*pi + 4*pi*pi + pi*pi/4 + pi*pi/9);
	
	EXPECT_MPREAL_EQ(pi/norme, quad.ComponantReal());
	EXPECT_MPREAL_EQ(2*pi/norme, quad.ComponantI());
	EXPECT_MPREAL_EQ(pi/2/norme, quad.ComponantJ());
	EXPECT_MPREAL_EQ(pi/3/norme, quad.ComponantK());
}

TEST_F(QuaternionTest,FormatterParser){
	Quaternion a(pi,2*pi,pi/2,pi/3);
	a.Normalize();
	stringstream sstr;
	sstr << a;
	auto text = sstr.str();
	std::stringstream sstr2(text);
	Quaternion quad;
	sstr2 >> quad;
	
	EXPECT_MPREAL_EQ(a.ComponantReal(), quad.ComponantReal());
	EXPECT_MPREAL_EQ(a.ComponantI(), quad.ComponantI());
	EXPECT_MPREAL_EQ(a.ComponantJ(), quad.ComponantJ());
	EXPECT_MPREAL_EQ(a.ComponantK(), quad.ComponantK());
}




