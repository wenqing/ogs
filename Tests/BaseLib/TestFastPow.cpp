/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestFastPow.cpp
 *  Created on July 20, 2018, 10:35 AM
 */

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <cstdlib>
#include "Tests/TestTools.h"

#include "BaseLib/Algorithm.h"

TEST(BaseLib, checkFastPow)
{
    const double min_a = std::numeric_limits<double>::epsilon();
    const double max_a = 1.0 / min_a;
    const double a =
        ((double(rand()) / double(RAND_MAX)) * (max_a - min_a)) + min_a;
    // Generate a random real number in [-16.0, 16.0]
    const double n = ((double(rand()) / double(RAND_MAX)) * 32.0) - 16.0;
    ASSERT_NEAR(std::pow(a, n), BaseLib::fastPow(a, n), min_a);
}
