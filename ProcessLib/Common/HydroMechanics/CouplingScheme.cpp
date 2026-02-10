// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CouplingScheme.h"

namespace ProcessLib
{
namespace Common::HydroMechanics
{
CouplingScheme parseCouplingScheme(
    std::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return Monolithic{};
    }

    auto const coupling_scheme_type =
        //! \ogs_file_param{prj__processes__process__common__HydroMechanics__coupling_scheme__type}
        config->getConfigParameter<std::string>("type");

    if (coupling_scheme_type == "monolithic")
    {
        return Monolithic{};
    }

    // Default value is 0.5, which is recommended by [Mikelic & Wheeler].
    double const fixed_stress_stabilization_parameter =
        //! \ogs_file_param{prj__processes__process__common__HydroMechanics__coupling_scheme__fixed_stress_stabilization_parameter}
        config->getConfigParameter<double>(
            "fixed_stress_stabilization_parameter", 0.5);

    DBUG("Using value {:g} for coupling parameter of staggered scheme.",
         fixed_stress_stabilization_parameter);

    {  // Check parameter value. Optimum is not a-priori known, but within
       // certain interval [Storvik & Nordbotten]
        double const csp_min = 1.0 / 6.0;
        double const csp_max = 1.0;
        if (fixed_stress_stabilization_parameter < csp_min ||
            fixed_stress_stabilization_parameter > csp_max)
        {
            WARN(
                "Value of coupling scheme parameter = {:g} is out of "
                "reasonable range ({:g}, {:g}).",
                fixed_stress_stabilization_parameter, csp_min, csp_max);
        }
    }

    bool const fixed_stress_over_time_step =
        //! \ogs_file_param{prj__processes__process__common__HydroMechanics__coupling_scheme__fixed_stress_over_time_step}
        config->getConfigParameter<std::string>("fixed_stress_over_time_step",
                                                "false") == "true";

    return Staggered{fixed_stress_stabilization_parameter,
                     fixed_stress_over_time_step};
}

}  // namespace Common::HydroMechanics
}  // namespace ProcessLib