#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include "Precision.h"
#include <Solid.h>
#include <SolidFormatter.h>
#include <SolidParser.h>

using namespace GeometricalSpaceObjects;
using namespace GeometricalSolid;

class AnyShape : public GeometricalSolid::Shape{
public:
	
	AnyShape() {
		invertedInertia = inertia.MatrixInverse();
	}
	
	virtual const double& Mass() const {
		return this->mass;
	}
	
	virtual const GeometricalSpaceObjects::Matrix<double>& InvertedIntertia() const {
		return this->invertedInertia;
	}

private:
	Matrix<double> inertia{pi,0,0,0,pi,0,0,0,pi};
	Matrix<double> invertedInertia;
	double mass{10.};
};

class SolidTest : public ::testing::Test {
public:
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid *s{nullptr};
protected:
	virtual void SetUp() {
		s = new Solid(std::move(shape));
		double dt = 0.0001;
		Matrix<double> I(pi,0,0,
										 0,pi,0,
										 0,0,pi);
		
		Vector<double> f(pi,pi/2,pi/4);
		Vector<double> M(pi,-pi/2,pi/5);
		s->Force(f);
		s->Momentum(M);
		s->UpdateVelocities(dt);
		s->UpdatePosition(dt);
	}
	virtual void TearDown() {}
};


TEST_F(SolidTest,DefaultConstructor){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	
	EXPECT_MPREAL_EQ(0, solid.Basis().Origin().CoordinateX());
	EXPECT_MPREAL_EQ(0, solid.Basis().Origin().CoordinateY());
	EXPECT_MPREAL_EQ(0, solid.Basis().Origin().CoordinateZ());
	
	EXPECT_MPREAL_EQ(1, solid.Basis().Orientation().ComponantReal());
	EXPECT_MPREAL_EQ(0, solid.Basis().Orientation().ComponantI());
	EXPECT_MPREAL_EQ(0, solid.Basis().Orientation().ComponantJ());
	EXPECT_MPREAL_EQ(0, solid.Basis().Orientation().ComponantK());
	
	EXPECT_MPREAL_EQ(0, solid.Velocity().ComponantX());
	EXPECT_MPREAL_EQ(0, solid.Velocity().ComponantY());
	EXPECT_MPREAL_EQ(0, solid.Velocity().ComponantZ());
	
	EXPECT_MPREAL_EQ(0, solid.AngularVelocity().ComponantX());
	EXPECT_MPREAL_EQ(0, solid.AngularVelocity().ComponantY());
	EXPECT_MPREAL_EQ(0, solid.AngularVelocity().ComponantZ());
	
	EXPECT_MPREAL_EQ(0, solid.Force().ComponantX());
	EXPECT_MPREAL_EQ(0, solid.Force().ComponantY());
	EXPECT_MPREAL_EQ(0, solid.Force().ComponantZ());
	
	EXPECT_MPREAL_EQ(0, solid.Momentum().ComponantX());
	EXPECT_MPREAL_EQ(0, solid.Momentum().ComponantY());
	EXPECT_MPREAL_EQ(0, solid.Momentum().ComponantZ());
	
	EXPECT_FALSE(solid.IsXTranslationLocked());
	EXPECT_FALSE(solid.IsYTranslationLocked());
	EXPECT_FALSE(solid.IsZTranslationLocked());
	
	EXPECT_FALSE(solid.IsXRotationLocked());
	EXPECT_FALSE(solid.IsYRotationLocked());
	EXPECT_FALSE(solid.IsZRotationLocked());
	
	EXPECT_MPREAL_EQ(10, solid.Shape()->Mass());
	
	auto m = solid.Shape()->InvertedIntertia();
	for(auto i = 0 ; i < 3 ; ++i) {
		for(auto j = 0 ; j < 3 ; ++j) {
			if(i != j)
				EXPECT_MPREAL_EQ(0, m.Element(i, j));
			else
				EXPECT_MPREAL_EQ(1./M_PI, m.Element(i, j));
		}
	}
}


TEST_F(SolidTest,UpdatePositionTranslation){
	
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));

	Vector<double> v(pi,pi/2,pi/4);
	double dt = 0.0001;
	solid.Velocity(v);
	solid.UpdatePosition(dt);
	
	EXPECT_MPREAL_EQ(solid.Basis().Origin().CoordinateX() , pi  *dt);
	EXPECT_MPREAL_EQ(solid.Basis().Origin().CoordinateY() , pi/2*dt);
	EXPECT_MPREAL_EQ(solid.Basis().Origin().CoordinateZ() , pi/4*dt);
}

TEST_F(SolidTest,UpdateVelocityTranslation){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));

	Vector<double> f(pi,pi/2,pi/4);
	double mass = 10.;
	double dt = 0.0001;
	solid.Force(f);
	solid.UpdateVelocities(dt);
	
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantX() , pi  *dt/mass);
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantY() , pi/2*dt/mass);
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantZ() , pi/4*dt/mass);
}

TEST_F(SolidTest,UpdateVelocityAndPositionTranslation){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> f(pi,pi/2,pi/4);
	double mass = 10.;
	double dt = 0.0001;
	solid.Force(f);
	
	solid.UpdateVelocities(dt);
	solid.UpdatePosition(dt);
	
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantX() , pi  *dt/mass);
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantY() , pi/2*dt/mass);
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantZ() , pi/4*dt/mass);
	EXPECT_MPREAL_EQ( solid.Basis().Origin().CoordinateX() , pi  *dt*dt/mass);
	EXPECT_MPREAL_EQ( solid.Basis().Origin().CoordinateY() , pi/2*dt*dt/mass);
	EXPECT_MPREAL_EQ( solid.Basis().Origin().CoordinateZ() , pi/4*dt*dt/mass);
}


TEST_F(SolidTest,UpdateVelocityAndPositionTranslationThetaMethod){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> f(pi,pi/2,pi/4);
	double mass = 10.;
	double dt = 0.0001;
	double T = 0;
	solid.Force(f);
	
	for(int i = 0 ; i < 10 ; i++){
		T += dt;
		solid.UpdatePosition(dt/2.);
		solid.UpdateVelocities(dt);
		solid.UpdatePosition(dt/2.);
	}
	
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantX() , pi  *T/mass);
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantY() , pi/2*T/mass);
	EXPECT_MPREAL_EQ( solid.Velocity().ComponantZ() , pi/4*T/mass);
	EXPECT_MPREAL_EQ( solid.Basis().Origin().CoordinateX(), pi  *T*T/2./mass);
	EXPECT_MPREAL_EQ( solid.Basis().Origin().CoordinateY() , pi/2*T*T/2./mass);
	EXPECT_MPREAL_EQ( solid.Basis().Origin().CoordinateZ() , pi/4*T*T/2./mass);
}


TEST_F(SolidTest,UpdatePositionRotationX){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> w(pi/6,0.,0.);
	double dt = 1;
	
	solid.AngularVelocity(w);
	solid.UpdatePosition(dt);
	
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantX() , pi/6);
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantReal() , cos(pi/12));
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantI() , sin(pi/12));
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantJ() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantK() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantX() , 1);
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantY() ,cos(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantZ() ,sin(pi/6));
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantY() ,-sin(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantZ() , cos(pi/6));
}


TEST_F(SolidTest,UpdatePositionRotationY){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> w(0 ,pi/6, 0);
	double dt = 1;
	
	solid.AngularVelocity(w);
	solid.UpdatePosition(dt);
	
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantY() , pi/6);
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantReal() , cos(pi/12));
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantI() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantJ() , sin(pi/12));
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantK() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantX() , cos(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantZ() ,-sin(pi/6));
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantY() , 1);
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantX() , sin(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantZ() , cos(pi/6));
}

TEST_F(SolidTest,UpdatePositionRotationZ){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> w(0 , 0, pi/6);
	double dt = 1;
	
	solid.AngularVelocity(w);
	solid.UpdatePosition(dt);
	
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.AngularVelocity().ComponantZ() , pi/6);
	
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantReal() , cos(pi/12));
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantI() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantJ() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().Orientation().ComponantK() , sin(pi/12));
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantX() , cos(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantY() , sin(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantX() ,-sin(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantY() , cos(pi/6));
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantZ() ,1);
}

TEST_F(SolidTest,UpdatePositionRotationXYZ){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> w;
	double dt = 1;
	
	// Rotation X
	w.SetComponants(pi,0,0);
	solid.AngularVelocity(w);
	solid.UpdatePosition(dt);
	// Rotation Y
	w.SetComponants(0,pi,0);
	solid.AngularVelocity(w);
	solid.UpdatePosition(dt);
	// Rotation Z
	w.SetComponants(0,0,pi);
	solid.AngularVelocity(w);
	solid.UpdatePosition(dt);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantX() , 1);
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisX().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantY() , 1);
	EXPECT_MPREAL_EQ( solid.Basis().AxisY().ComponantZ() , 0);
	
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantX() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantY() , 0);
	EXPECT_MPREAL_EQ( solid.Basis().AxisZ().ComponantZ() , 1);
}

TEST_F(SolidTest,UpdateAngularVelocity){
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	Vector<double> M(pi,-pi/2,pi/5);
	Matrix<double> I(pi,0,0,
					 0,pi,0,
					 0,0,pi);
	
	double dt = 1;
	solid.Momentum(M);
	solid.UpdateVelocities(dt);
	
	EXPECT_MPREAL_EQ(solid.AngularVelocity().ComponantX(),M.ComponantX()*dt/pi);
	EXPECT_MPREAL_EQ(solid.AngularVelocity().ComponantY(),M.ComponantY()*dt/pi);
	EXPECT_MPREAL_EQ(solid.AngularVelocity().ComponantZ(),M.ComponantZ()*dt/pi);
}


TEST_F(SolidTest,Formatter){
	std::string expected ("3.1415926535897932e-09\t1.5707963267948966e-09\t7.8539816339744829e-10\n1.0000000000000000e+00\t5.0000000000000001e-09\t-2.5000000000000001e-09	1.0000000000000001e-09\n3.1415926535897928e-05\t1.5707963267948964e-05\t7.8539816339744820e-06\n1.0000000000000000e-04\t-5.0000000000000002e-05\t2.0000000000000002e-05\n3.1415926535897931e+00\t1.5707963267948966e+00\t7.8539816339744828e-01\n3.1415926535897931e+00\t-1.5707963267948966e+00\t6.2831853071795862e-01");
	std::stringstream sstr;
	sstr << *s;
	auto text = sstr.str();
	EXPECT_TRUE(expected == text);
}

TEST_F(SolidTest,Parser){
	std::string toParse ("3.1415926535897932e-09\t1.5707963267948966e-09\t7.8539816339744829e-10\n1.0000000000000000e+00\t5.0000000000000001e-09\t-2.5000000000000001e-09	1.0000000000000001e-09\n3.1415926535897928e-05\t1.5707963267948964e-05\t7.8539816339744820e-06\n1.0000000000000000e-04\t-5.0000000000000002e-05\t2.0000000000000002e-05\n3.1415926535897931e+00\t1.5707963267948966e+00\t7.8539816339744828e-01\n3.1415926535897931e+00\t-1.5707963267948966e+00\t6.2831853071795862e-01");
	std::stringstream sstr(toParse);
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	sstr >> solid;
	
	EXPECT_MPREAL_EQ(s->Basis().Origin().CoordinateX(), solid.Basis().Origin().CoordinateX());
	EXPECT_MPREAL_EQ(s->Basis().Origin().CoordinateY(), solid.Basis().Origin().CoordinateY());
	EXPECT_MPREAL_EQ(s->Basis().Origin().CoordinateZ(), solid.Basis().Origin().CoordinateZ());
	
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantReal(), solid.Basis().Orientation().ComponantReal());
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantI(), solid.Basis().Orientation().ComponantI());
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantJ(), solid.Basis().Orientation().ComponantJ());
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantK(), solid.Basis().Orientation().ComponantK());
	
	EXPECT_MPREAL_EQ(s->Velocity().ComponantX(), solid.Velocity().ComponantX());
	EXPECT_MPREAL_EQ(s->Velocity().ComponantY(), solid.Velocity().ComponantY());
	EXPECT_MPREAL_EQ(s->Velocity().ComponantZ(), solid.Velocity().ComponantZ());
	
	EXPECT_MPREAL_EQ(s->AngularVelocity().ComponantX(), solid.AngularVelocity().ComponantX());
	EXPECT_MPREAL_EQ(s->AngularVelocity().ComponantY(), solid.AngularVelocity().ComponantY());
	EXPECT_MPREAL_EQ(s->AngularVelocity().ComponantZ(), solid.AngularVelocity().ComponantZ());
	
	EXPECT_MPREAL_EQ(s->Force().ComponantX(), solid.Force().ComponantX());
	EXPECT_MPREAL_EQ(s->Force().ComponantY(), solid.Force().ComponantY());
	EXPECT_MPREAL_EQ(s->Force().ComponantZ(), solid.Force().ComponantZ());
	
	EXPECT_MPREAL_EQ(s->Momentum().ComponantX(), solid.Momentum().ComponantX());
	EXPECT_MPREAL_EQ(s->Momentum().ComponantY(), solid.Momentum().ComponantY());
	EXPECT_MPREAL_EQ(s->Momentum().ComponantZ(), solid.Momentum().ComponantZ());
}

TEST_F(SolidTest,FormatterParser){
	std::stringstream sstr;
	sstr << *s;
	auto text = sstr.str();

	std::stringstream sstr2(text);
	std::unique_ptr<Shape> shape{new AnyShape()};
	Solid solid(std::move(shape));
	sstr2 >> solid;
	
	EXPECT_MPREAL_EQ(s->Basis().Origin().CoordinateX(), solid.Basis().Origin().CoordinateX());
	EXPECT_MPREAL_EQ(s->Basis().Origin().CoordinateY(), solid.Basis().Origin().CoordinateY());
	EXPECT_MPREAL_EQ(s->Basis().Origin().CoordinateZ(), solid.Basis().Origin().CoordinateZ());
	
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantReal(), solid.Basis().Orientation().ComponantReal());
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantI(), solid.Basis().Orientation().ComponantI());
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantJ(), solid.Basis().Orientation().ComponantJ());
	EXPECT_MPREAL_EQ(s->Basis().Orientation().ComponantK(), solid.Basis().Orientation().ComponantK());
	
	EXPECT_MPREAL_EQ(s->Velocity().ComponantX(), solid.Velocity().ComponantX());
	EXPECT_MPREAL_EQ(s->Velocity().ComponantY(), solid.Velocity().ComponantY());
	EXPECT_MPREAL_EQ(s->Velocity().ComponantZ(), solid.Velocity().ComponantZ());
	
	EXPECT_MPREAL_EQ(s->AngularVelocity().ComponantX(), solid.AngularVelocity().ComponantX());
	EXPECT_MPREAL_EQ(s->AngularVelocity().ComponantY(), solid.AngularVelocity().ComponantY());
	EXPECT_MPREAL_EQ(s->AngularVelocity().ComponantZ(), solid.AngularVelocity().ComponantZ());
	
	EXPECT_MPREAL_EQ(s->Force().ComponantX(), solid.Force().ComponantX());
	EXPECT_MPREAL_EQ(s->Force().ComponantY(), solid.Force().ComponantY());
	EXPECT_MPREAL_EQ(s->Force().ComponantZ(), solid.Force().ComponantZ());
	
	
	EXPECT_MPREAL_EQ(s->Momentum().ComponantX(), solid.Momentum().ComponantX());
	EXPECT_MPREAL_EQ(s->Momentum().ComponantY(), solid.Momentum().ComponantY());
	EXPECT_MPREAL_EQ(s->Momentum().ComponantZ(), solid.Momentum().ComponantZ());
}
