/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/CurveProperty.h"

namespace MaterialPropertyLib
{
CurveProperty::CurveProperty(Variable const type,
                             MathLib::PiecewiseLinearInterpolation const& curve)
    : _variable_type(type), _curve(curve)
{
}

PropertyDataType CurveProperty::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    auto const x =
        std::get<double>(variable_array[static_cast<int>(_variable_type)]);
    return _curve.getValue(x);
}

PropertyDataType CurveProperty::dValue(
    VariableArray const& variable_array,
    Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    auto const x =
        std::get<double>(variable_array[static_cast<int>(_variable_type)]);
    return _curve.getDerivative(x);
}
}  // namespace MaterialPropertyLib
