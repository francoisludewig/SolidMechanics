#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

#include <Basis.h>
#include <Point.h>
#include "Precision.h"


using namespace std;
using namespace GeometricalSpaceObjects;


class BasisTest : public ::testing::Test {
public:
	Basis<Type> basis;

protected:
	virtual void SetUp() {
#ifndef DOUBLE_PRECISON
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(50));
#endif
		
		basis.Origin(Point<Type>(pi, pi/4, 5*pi));
		Quaternion<Type> quad(pi, pi/7, 6*pi, pi/11);
		quad.Normalize();
		basis.Orientation(std::move(quad));
	}
	
	virtual void TearDown() {}
};



TEST_F(BasisTest,Constructor){
	Basis<Type> a;
	Quaternion<Type> q = a.Orientation();
	
	EXPECT_MPREAL_EQ(1,q.ComponantReal());
	EXPECT_MPREAL_EQ(0,q.ComponantI());
	EXPECT_MPREAL_EQ(0,q.ComponantJ());
	EXPECT_MPREAL_EQ(0,q.ComponantK());
	
	Point<Type> o = a.Origin();
	EXPECT_MPREAL_EQ(0,o.CoordinateX());
	EXPECT_MPREAL_EQ(0,o.CoordinateY());
	EXPECT_MPREAL_EQ(0,o.CoordinateZ());
	
	
	Vector<Type> e1 = a.AxisX();
	EXPECT_MPREAL_EQ(1,e1.ComponantX());
	EXPECT_MPREAL_EQ(0,e1.ComponantY());
	EXPECT_MPREAL_EQ(0,e1.ComponantZ());
	
	Vector<Type> e2 = a.AxisY();
	EXPECT_MPREAL_EQ(0,e2.ComponantX());
	EXPECT_MPREAL_EQ(1,e2.ComponantY());
	EXPECT_MPREAL_EQ(0,e2.ComponantZ());
	
	Vector<Type> e3 = a.AxisZ();
	EXPECT_MPREAL_EQ(0,e3.ComponantX());
	EXPECT_MPREAL_EQ(0,e3.ComponantY());
	EXPECT_MPREAL_EQ(1,e3.ComponantZ());
}

TEST_F(BasisTest,translate){
	Basis<Type> a;
	Vector<Type> o;
	o.SetComponants(pi,-pi/3,7*pi/6);
	a += o;
	
	EXPECT_MPREAL_EQ(a.Origin().CoordinateX() , o.ComponantX() );
	EXPECT_MPREAL_EQ(a.Origin().CoordinateY() , o.ComponantY() );
	EXPECT_MPREAL_EQ(a.Origin().CoordinateZ() , o.ComponantZ() );
}


TEST_F(BasisTest,rotate){
	Basis<Type> a,b;
	Vector<Type> w;
	w.SetComponants(0,pi,0);
	Quaternion<Type> q(w);
	
	b = a*q;
	
	EXPECT_MPREAL_EQ(b.AxisX().ComponantX() , -a.AxisX().ComponantX() );
	EXPECT_MPREAL_EQ(b.AxisX().ComponantY() , -a.AxisX().ComponantY() );
	EXPECT_MPREAL_EQ(b.AxisX().ComponantZ() , -a.AxisX().ComponantZ() );
	EXPECT_MPREAL_EQ(b.AxisY().ComponantX() , a.AxisY().ComponantX() );
	EXPECT_MPREAL_EQ(b.AxisY().ComponantY() , a.AxisY().ComponantY() );
	EXPECT_MPREAL_EQ(b.AxisY().ComponantZ() , a.AxisY().ComponantZ() );
	EXPECT_MPREAL_EQ(b.AxisZ().ComponantX() , -a.AxisZ().ComponantX() );
	EXPECT_MPREAL_EQ(b.AxisZ().ComponantY() , -a.AxisZ().ComponantY() );
	EXPECT_MPREAL_EQ(b.AxisZ().ComponantZ() , -a.AxisZ().ComponantZ() );
	
	a*=q;
	
	EXPECT_MPREAL_EQ(b.AxisX().ComponantX() , a.AxisX().ComponantX() );
	EXPECT_MPREAL_EQ(b.AxisX().ComponantY() , a.AxisX().ComponantY() );
	EXPECT_MPREAL_EQ(b.AxisX().ComponantZ() , a.AxisX().ComponantZ() );
	EXPECT_MPREAL_EQ(b.AxisY().ComponantX() , a.AxisY().ComponantX() );
	EXPECT_MPREAL_EQ(b.AxisY().ComponantY() , a.AxisY().ComponantY() );
	EXPECT_MPREAL_EQ(b.AxisY().ComponantZ() , a.AxisY().ComponantZ() );
	EXPECT_MPREAL_EQ(b.AxisZ().ComponantX() , a.AxisZ().ComponantX() );
	EXPECT_MPREAL_EQ(b.AxisZ().ComponantY() , a.AxisZ().ComponantY() );
	EXPECT_MPREAL_EQ(b.AxisZ().ComponantZ() , a.AxisZ().ComponantZ() );
}


TEST_F(BasisTest,Construct){
	Basis<Type> a;
	Vector<Type> u(2*pi,pi/2,pi/3);
	a.AxisX(u);
	Vector<Type> n,t,s;
	n = a.AxisX();
	t = a.AxisY();
	s = a.AxisZ();
	
	EXPECT_MPREAL_EQ(0,n*t);
	EXPECT_MPREAL_EQ(0,n*s);
	EXPECT_MPREAL_EQ(0,s*t);
	
	EXPECT_MPREAL_EQ(1,n.Norme());
	EXPECT_MPREAL_EQ(1,t.Norme());
	EXPECT_MPREAL_EQ(1,s.Norme());
	
}


TEST_F(BasisTest,Formatter){
	auto text = LuGaBasisFormatter<Type>().Format(basis);
	
	std::cout << text << std::endl;
	
	EXPECT_TRUE(std::string("3.141592653589793e+00	7.853981633974483e-01	1.570796326794897e+01\n1.643353249699871e-01	2.347647499571244e-02	9.860119498199226e-01	1.493957499727155e-02") == text);
}

TEST_F(BasisTest,Parser){
	auto ba = LuGaBasisParser<Type>().Parse(std::string("3.141592653589793e+00	7.853981633974483e-01	1.570796326794897e+01\n1.643353249699871e-01	2.347647499571244e-02	9.860119498199226e-01	1.493957499727155e-02"));
	
	EXPECT_MPREAL_EQ(basis.Origin().CoordinateX(), ba.Origin().CoordinateX());
	EXPECT_MPREAL_EQ(basis.Origin().CoordinateY(), ba.Origin().CoordinateY());
	EXPECT_MPREAL_EQ(basis.Origin().CoordinateZ(), ba.Origin().CoordinateZ());
	
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantReal(), ba.Orientation().ComponantReal());
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantI(), ba.Orientation().ComponantI());
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantJ(), ba.Orientation().ComponantJ());
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantK(), ba.Orientation().ComponantK());
}

TEST_F(BasisTest,FormatterParser){
	auto text = LuGaBasisFormatter<Type>().Format(basis);
	auto ba = LuGaBasisParser<Type>().Parse(text);
	EXPECT_MPREAL_EQ(basis.Origin().CoordinateX(), ba.Origin().CoordinateX());
	EXPECT_MPREAL_EQ(basis.Origin().CoordinateY(), ba.Origin().CoordinateY());
	EXPECT_MPREAL_EQ(basis.Origin().CoordinateZ(), ba.Origin().CoordinateZ());
	
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantReal(), ba.Orientation().ComponantReal());
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantI(), ba.Orientation().ComponantI());
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantJ(), ba.Orientation().ComponantJ());
	EXPECT_MPREAL_EQ(basis.Orientation().ComponantK(), ba.Orientation().ComponantK());
}


TEST_F(BasisTest,IO_Operator){
	Basis<Type> a;
	Basis<Type> b;
	Quaternion<Type> q(pi,2*pi,pi/2,pi/3);
	Point<Type> p(2*pi,pi/2,pi/3);
	a.Orientation(q);
	a.Origin(p);
	
	ofstream fichierOut("testBasis.txt", ios::out | ios::trunc);
	if(fichierOut){
		//		fichierOut << a;
		fichierOut.close();
	}
	
	ifstream fichierIn("testBasis.txt", ios::in);
	if(fichierIn){
		//fichierIn >> b;
		fichierIn.close();
	}
	
	Point<Type> o_a = a.Origin();
	Point<Type> o_b = b.Origin();
	
	EXPECT_MPREAL_EQ(o_a.CoordinateX(),o_b.CoordinateX());
	EXPECT_MPREAL_EQ(o_a.CoordinateY(),o_b.CoordinateY());
	EXPECT_MPREAL_EQ(o_a.CoordinateZ(),o_b.CoordinateZ());
	
	Quaternion<Type> q_a = a.Orientation();
	Quaternion<Type> q_b = b.Orientation();
	
	EXPECT_MPREAL_EQ(q_a.ComponantReal(),q_b.ComponantReal());
	EXPECT_MPREAL_EQ(q_a.ComponantI(),q_b.ComponantI());
	EXPECT_MPREAL_EQ(q_a.ComponantJ(),q_b.ComponantJ());
	EXPECT_MPREAL_EQ(q_a.ComponantK(),q_b.ComponantK());
	
	remove("testBasis.txt");
}
