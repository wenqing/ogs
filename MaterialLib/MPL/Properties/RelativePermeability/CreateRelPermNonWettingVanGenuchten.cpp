/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 14, 2020, 3:50 PM
 */

#include "CreateRelPermNonWettingVanGenuchten.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "RelPermNonWettingVanGenuchten.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createRelPermNonWettingVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "RelativePermeabilityNonWettingVanGenuchten");
    DBUG("Create RelPermNonWettingVanGenuchten medium property");

    auto const residual_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingVanGenuchten__residual_saturation}
        config.getConfigParameter<double>("residual_saturation");
    auto const maximum_saturation =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingVanGenuchten__maximum_saturation}
        config.getConfigParameter<double>("maximum_saturation");

    auto const exponent =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");

    auto const min_relative_permeability =
        //! \ogs_file_param{properties__property__RelativePermeabilityNonWettingVanGenuchten__minimum_relative_permeability}
        config.getConfigParameter<double>("minimum_relative_permeability");

    return std::make_unique<RelPermNonWettingVanGenuchten>(
        residual_saturation, maximum_saturation, exponent,
        min_relative_permeability);
}
}  // namespace MaterialPropertyLib