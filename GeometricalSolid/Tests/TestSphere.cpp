#include <gtest/gtest.h>
#include "Precision.h"
#include <memory.h>
#include "Sphere.h"
#include "Solid.h"

using namespace GeometricalSolid;

class SphereTest : public ::testing::Test {
public:
	Sphere s;
protected:
	virtual void SetUp() {
		s = Sphere(0.01, 2500);
	}
	virtual void TearDown() {}
};

TEST_F(SphereTest,DefaultConstructor) {
	EXPECT_MPREAL_EQ(0.01, s.Radius());
	EXPECT_MPREAL_EQ(2500, s.Density());
	EXPECT_MPREAL_EQ(4./3.*M_PI*0.000001, s.Volume());
	EXPECT_MPREAL_EQ(2500.*4./3.*M_PI*0.000001, s.Mass());
	auto m = s.InvertedIntertia();
	for(auto i = 0 ; i < 3 ; ++i) {
		for(auto j = 0 ; j < 3 ; ++j) {
			if(i != j)
				EXPECT_MPREAL_EQ(0, m.Element(i, j));
			else
				EXPECT_MPREAL_EQ(5./(2*0.0001*2500.*4./3.*M_PI*0.000001), m.Element(i, j));
		}
	}
	
	EXPECT_TRUE(Shape::Nature::Particle == s.Nature());
	EXPECT_TRUE(Shape::Form::Sphere == s.Form());
}

TEST(SphereSolidTest,inSolid) {
	std::unique_ptr<Shape> sph(new Sphere(0.01, 2500));
	Solid solid(std::move(sph));
	Sphere sphere(0.01, 2500);
	const Sphere* sph2 = static_cast<const Sphere*>(solid.Shape());
	
	EXPECT_MPREAL_EQ(sph2->Radius(), sphere.Radius());
	EXPECT_MPREAL_EQ(sph2->Density(), sphere.Density());
	EXPECT_MPREAL_EQ(sph2->Volume(), sphere.Volume());
	EXPECT_MPREAL_EQ(sph2->Mass(), sphere.Mass());
	auto m = sphere.InvertedIntertia();
	auto m2 = sph2->InvertedIntertia();

	for(auto i = 0 ; i < 3 ; ++i) {
		for(auto j = 0 ; j < 3 ; ++j) {
				EXPECT_MPREAL_EQ(m2.Element(i, j), m.Element(i, j));
		}
	}
	
	EXPECT_TRUE(sph2->Nature() == sphere.Nature());
	EXPECT_TRUE(sph2->Form() == sphere.Form());
}
