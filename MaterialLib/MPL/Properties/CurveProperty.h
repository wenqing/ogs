/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialPropertyLib
{
/// The CurveProperty class. The property reads the value from a curve.
/// The current implementation accepts only the double datatype defined in
/// PropertyDataType.
class CurveProperty final : public Property
{
public:
    /// This constructor accepts a curve.
    CurveProperty(Variable const type,
                  MathLib::PiecewiseLinearInterpolation const& curve);

    /// This method computes the value of a property depending linearly on
    /// the value of the given primary variable.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// This method will compute the derivative of a property with respect to
    /// the given primary variable.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

private:
    Variable const _variable_type;
    MathLib::PiecewiseLinearInterpolation const& _curve;
};
}  // namespace MaterialPropertyLib
