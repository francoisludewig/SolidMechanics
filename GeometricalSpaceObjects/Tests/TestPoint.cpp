#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "Point.h"
#include "Vector.h"
#include "Precision.h"

using namespace std;
using namespace GeometricalSpaceObjects;


#define Vector Vector<Type>
#define Point Point<Type>


class PointTest : public ::testing::Test {
protected:
	virtual void SetUp() {		
#ifndef DOUBLE_PRECISON
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(50));
#endif
	}
	
	virtual void TearDown() {}
};



TEST_F(PointTest,Constructor){
	Point a;
	
	EXPECT_MPREAL_EQ(0,a.CoordinateX());
	EXPECT_MPREAL_EQ(0,a.CoordinateY());
	EXPECT_MPREAL_EQ(0,a.CoordinateZ());
}

TEST_F(PointTest,Constructor2){
	Point a(pi,2*pi,pi/5);
	
	EXPECT_MPREAL_EQ(pi, a.CoordinateX());
	EXPECT_MPREAL_EQ(2*pi,a.CoordinateY());
	EXPECT_MPREAL_EQ(pi/5,a.CoordinateZ());
}


TEST_F(PointTest,Calculus){
	Point a(pi   ,2*pi,pi/5);
	Point b(-3*pi,pi  ,4*pi/5);
	
	Vector c = b-a;
	
	EXPECT_MPREAL_EQ(-4*pi , c.ComponantX());
	EXPECT_MPREAL_EQ(-pi   , c.ComponantY());
	EXPECT_MPREAL_EQ(3*pi/5, c.ComponantZ());
	
	Point d = a + c;
	
	EXPECT_MPREAL_EQ(d.CoordinateX(), b.CoordinateX());
	EXPECT_MPREAL_EQ(d.CoordinateY(), b.CoordinateY());
	EXPECT_MPREAL_EQ(d.CoordinateZ(), b.CoordinateZ());
	
	a += c;
	
	EXPECT_MPREAL_EQ(a.CoordinateX(), b.CoordinateX());
	EXPECT_MPREAL_EQ(a.CoordinateY(), b.CoordinateY());
	EXPECT_MPREAL_EQ(a.CoordinateZ(), b.CoordinateZ());
}


TEST_F(PointTest,IO_Operator){
	Point a(pi,2*pi,pi/5);
	Point b;
	
	ofstream fichierOut("testPoint.txt", ios::out | ios::trunc);
	if(fichierOut){
		fichierOut << a;
		fichierOut.close();
	}
	
	ifstream fichierIn("testPoint.txt", ios::in);
	if(fichierIn){
		fichierIn >> b;
		fichierIn.close();
	}
	
	EXPECT_MPREAL_EQ(b.CoordinateX(), a.CoordinateX());
	EXPECT_MPREAL_EQ(b.CoordinateY(), a.CoordinateY());
	EXPECT_MPREAL_EQ(b.CoordinateZ(), a.CoordinateZ());
	
	remove("testPoint.txt");
}

TEST_F(PointTest,Formatter){
	Point a(pi,2*pi,pi/5);
	std::stringstream sstr;
	sstr << a;
	auto text = sstr.str();
	EXPECT_TRUE(std::string("3.1415926535897931e+00	6.2831853071795862e+00	6.2831853071795862e-01") == text);
}


TEST_F(PointTest,Parser){
	std::stringstream sstr(std::string("3.1415926535897931e+00	6.2831853071795862e+00	6.2831853071795862e-01"));
	Point pt;
	sstr >> pt;
	EXPECT_MPREAL_EQ(pi, pt.CoordinateX());
	EXPECT_MPREAL_EQ(2*pi, pt.CoordinateY());
	EXPECT_MPREAL_EQ(pi/5, pt.CoordinateZ());
}

TEST_F(PointTest,FormatterParser){
	Point a(pi,2*pi,pi/2);
	std::stringstream sstr;
	sstr << a;
	auto text = sstr.str();
	std::stringstream sstr2(text);
	Point pt;
	sstr2 >> pt;
	EXPECT_MPREAL_EQ(a.CoordinateX(), pt.CoordinateX());
	EXPECT_MPREAL_EQ(a.CoordinateY(), pt.CoordinateY());
	EXPECT_MPREAL_EQ(a.CoordinateZ(), pt.CoordinateZ());
}




