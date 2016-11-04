/**
 * \file
 * \author Thomas Fischer
 * \date   2010-09-07
 * \brief  Implementation of the PiecewiseLinearInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

#include "PiecewiseLinearInterpolation.h"

#include "BaseLib/quicksort.h"

namespace MathLib
{
PiecewiseLinearInterpolation::PiecewiseLinearInterpolation(
    std::vector<double>&& supporting_points,
    std::vector<double>&& values_at_supp_pnts,
    bool supp_pnts_sorted)
    : _supp_pnts(std::move(supporting_points)),
      _values_at_supp_pnts(std::move(values_at_supp_pnts))
{
    if (!supp_pnts_sorted)
    {
        BaseLib::quicksort<double, double>(
            _supp_pnts, static_cast<std::size_t>(0), _supp_pnts.size(),
            _values_at_supp_pnts);
    }
}

double PiecewiseLinearInterpolation::interpolate(
    std::vector<double> const& supp_pnts,
    std::vector<double> const& values_at_supp_pnts,
    std::size_t const interval_idx,
    double const pnt_to_interpolate) const
{
    assert(std::fabs(supp_pnts[interval_idx + 1] - supp_pnts[interval_idx]) >
           std::numeric_limits<double>::epsilon());

    // compute linear interpolation polynom: y = m * (x - support[i]) + value[i]
    const double m((values_at_supp_pnts[interval_idx + 1] -
                    values_at_supp_pnts[interval_idx]) /
                   (supp_pnts[interval_idx + 1] - supp_pnts[interval_idx]));

    return m * (pnt_to_interpolate - supp_pnts[interval_idx]) +
           values_at_supp_pnts[interval_idx];
}

double PiecewiseLinearInterpolation::getValue(
    double const pnt_to_interpolate) const
{
    if (pnt_to_interpolate <= _supp_pnts.front())
    {
        return _values_at_supp_pnts[0];
    }

    if (_supp_pnts.back() <= pnt_to_interpolate)
    {
        return _values_at_supp_pnts[_supp_pnts.size() - 1];
    }

    // search interval that has the point inside
    auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(),
                                    pnt_to_interpolate));
    std::size_t const interval_idx = std::distance(_supp_pnts.begin(), it) - 1;

    return interpolate(_supp_pnts, _values_at_supp_pnts, interval_idx,
                       pnt_to_interpolate);
}

double PiecewiseLinearInterpolation::getInverseValue(double const value) const
{
    std::size_t interval_idx = 0;
    if (_values_at_supp_pnts.front() < _values_at_supp_pnts.back())
    {
        if (value <= _values_at_supp_pnts.front())
        {
            return _supp_pnts[0];
        }
        else if (value >= _values_at_supp_pnts.back())
        {
            return _supp_pnts[_supp_pnts.size() - 1];
        }
        else
        {
            // search interval that has the point inside
            auto const& it(std::lower_bound(_values_at_supp_pnts.begin(),
                                            _values_at_supp_pnts.end(), value));
            interval_idx = std::distance(_values_at_supp_pnts.begin(), it) - 1;
        }
    }
    else
    {
        if (value >= _values_at_supp_pnts.front())
        {
            return _supp_pnts[0];
        }
        else if (value <= _values_at_supp_pnts.back())
        {
            return _supp_pnts[_supp_pnts.size() - 1];
        }
        else
        {
            // search interval in the reverse direction for the point inside
            auto const& it(std::lower_bound(_values_at_supp_pnts.rbegin(),
                                            _values_at_supp_pnts.rend(),
                                            value));
            interval_idx = _values_at_supp_pnts.size() -
                           (std::distance(_values_at_supp_pnts.rbegin(), it)) -
                           1;
        }
    }

    return interpolate(_values_at_supp_pnts, _supp_pnts, interval_idx, value);
}

double PiecewiseLinearInterpolation::getDerivative(
    double const pnt_to_interpolate) const
{
    if (pnt_to_interpolate < _supp_pnts.front() ||
        _supp_pnts.back() < pnt_to_interpolate)
    {
        return 0;
    }

    auto const& it(std::lower_bound(_supp_pnts.begin(), _supp_pnts.end(),
                                    pnt_to_interpolate));
    std::size_t interval_idx = std::distance(_supp_pnts.begin(), it);

    if (pnt_to_interpolate == _supp_pnts.front())
    {
        interval_idx = 1;
    }

    if (interval_idx > 2 && interval_idx < _supp_pnts.size() - 1)
    {
        double const tangent_right =
            (_values_at_supp_pnts[interval_idx - 1] -
             _values_at_supp_pnts[interval_idx + 1]) /
            (_supp_pnts[interval_idx - 1] - _supp_pnts[interval_idx + 1]);
        double const tangent_left =
            (_values_at_supp_pnts[interval_idx - 2] -
             _values_at_supp_pnts[interval_idx]) /
            (_supp_pnts[interval_idx - 2] - _supp_pnts[interval_idx]);
        double const w =
            (pnt_to_interpolate - _supp_pnts[interval_idx]) /
            (_supp_pnts[interval_idx - 1] - _supp_pnts[interval_idx]);
        return (1. - w) * tangent_right + w * tangent_left;
    }
    else
    {
        return (_values_at_supp_pnts[interval_idx] -
                _values_at_supp_pnts[interval_idx - 1]) /
               (_supp_pnts[interval_idx] - _supp_pnts[interval_idx - 1]);
    }
}

double PiecewiseLinearInterpolation::getSupportMax() const
{
    assert(!_supp_pnts.empty());
    return _supp_pnts.back();
}
double PiecewiseLinearInterpolation::getSupportMin() const
{
    assert(!_supp_pnts.empty());
    return _supp_pnts.front();
}

bool PiecewiseLinearInterpolation::isMonotonic() const
{
    const double tangent0 = getDerivative(_supp_pnts[0]);

    if (std::fabs(tangent0) < std::numeric_limits<double>::min())
        return false;
    else
    {
        return std::none_of(_supp_pnts.begin(), _supp_pnts.end(),
                            [&](const double p) {
                                return this->getDerivative(p) * tangent0 <= 0.;
                            });
    }
}

}  // end MathLib
