/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <algorithm>
#include <limits>
#include <range/v3/algorithm/min.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <vector>

#include "BaseLib/Error.h"
#include "cubic_roots.hpp"
namespace MathLib
{

class CubicSolver
{
public:
    CubicSolver(const double a, const double b, const double c, const double d)
        : a_(a), b_(b), c_(c), d_(d)
    {
        if (std::abs(a_) < 1e-9)
        {
            OGS_FATAL("'a' must be non-zero for a cubic equation.");
        }
    }

    // Method to solve the cubic equation using Boost
    std::vector<double> solve()
    {
        // Solve using Boost's cubic_roots
        std::array<double, 3> const roots =
            boost::math::tools::cubic_roots<double>(a_, b_, c_, d_);

        std::vector<double> filtered_roots =
            ranges::views::ref(roots) |
            ranges::views::filter([](double const root)
                                  { return !std::isnan(root); }) |
            ranges::to<std::vector>();
        ranges::sort(filtered_roots);

        return filtered_roots;
    }

    // Method that returns the smallest positive real root
    double smallestPositiveRealRoot()
    {
        std::vector<double> const roots = solve();

        auto positive_roots =
            roots | ranges::views::filter([](double root) { return root > 0; });

        // If no positive root exists, return NaN
        if (ranges::empty(positive_roots))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return ranges::min(positive_roots);
    }

private:
    double a_, b_, c_, d_;  // Coefficients of the cubic equation
};

}  // namespace MathLib