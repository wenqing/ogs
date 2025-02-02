/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveSetting.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::eval(
    ConstitutiveModels<DisplacementDim>& /*models*/, double const /*t*/,
    double const /*dt*/, ParameterLib::SpatialPosition const& /*x_position*/,
    MaterialPropertyLib::Medium const& /*medium*/,
    TemperatureData const /*T_data*/,
    CapillaryPressureData<DisplacementDim> const& /*p_cap_data*/,
    KelvinVector<DisplacementDim> const& /*eps_arg*/,
    StatefulData<DisplacementDim>& /*state*/,
    StatefulDataPrev<DisplacementDim> const& /*prev_state*/,
    ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>&
    /*mat_state*/,
    ConstitutiveTempData<DisplacementDim>& /*tmp*/,
    OutputData<DisplacementDim>& /*out*/,
    ConstitutiveData<DisplacementDim>& /*cd*/) const
{
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ProcessLib::RichardsMechanics
