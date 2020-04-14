/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#include "RelPermNonWettingVanGenuchten.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
RelPermNonWettingVanGenuchten::RelPermNonWettingVanGenuchten(
    const double Snr, const double Snmax, const double m, const double krel_min)
    : _saturation_r(1. - Snmax),
      _saturation_max(1. - Snr),
      _m(m),
      _krel_min(krel_min)
{
    if (!(_m > 0 && _m < 1))
    {
        OGS_FATAL(
            "The exponent value m = {:g} of van Genuchten relative "
            "permeability model, is out of its range of (0, 1)",
            _m);
    }
}

PropertyDataType RelPermNonWettingVanGenuchten::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_l = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);
    const double S = std::clamp(S_l, _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double krel =
        std::cbrt(1.0 - Se) * std::pow(1.0 - std::pow(Se, 1.0 / _m), 2.0 * _m);
    return std::max(_krel_min, krel);
}

PropertyDataType RelPermNonWettingVanGenuchten::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    const double S_l = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    const double S = std::clamp(S_l, _saturation_r + _minor_offset,
                                _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double cbrt1_Se = std::cbrt(1.0 - Se);
    const double temp_val = 1.0 - std::pow(Se, 1.0 / _m);
    return (-std::pow(temp_val, 2. * _m) / (3. * cbrt1_Se * cbrt1_Se) -
            2. * cbrt1_Se * std::pow(temp_val, 2. * _m - 1.) *
                std::pow(Se, (1. - _m) / _m)) /
           (_saturation_max - _saturation_r);
}

}  // namespace MaterialPropertyLib