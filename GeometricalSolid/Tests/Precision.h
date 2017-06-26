#pragma once

#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION

#define Type double
#define pi M_PI
#define EXPECT_MPREAL_EQ(x,y) EXPECT_TRUE(fabs(x-y) < std::numeric_limits<Type>::epsilon())

#else

#include "mpreal.h"
#define Type mpfr::mpreal
#define pi mpfr::const_pi()
#define EXPECT_MPREAL_EQ(x,y) EXPECT_TRUE(fabs(x-y) < std::numeric_limits<Type>::epsilon())

#endif

