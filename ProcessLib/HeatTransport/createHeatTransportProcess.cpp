/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createHeatTransportProcess.cpp
 *
 */
#include "createHeatTransportProcess.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "HeatTransportProcess.h"

namespace ProcessLib
{
namespace HeatTransport
{
std::unique_ptr<Process> createHeatTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "HEAT_TRANSPORT");

    DBUG("Create HeatTransportProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {//! \ogs_file_param_special{process__HEAT_TRANSPORT__process_variables__process_variable}
         "process_variable"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HeatTransport_temperature"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    //! \ogs_file_param{process__HEAT_TRANSPORT__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    // Parsing the solid properties
    //! \ogs_file_param{process__HEAT_TRANSPORT__material_property_solid}
    auto const& mat_solid_config = mat_config.getConfigSubtree("solid");

    // Density of solid.
    auto& solid_density = findParameter<double>(
        mat_solid_config,
        //! \ogs_file_param_special{process__HEAT_TRANSPORT__material_property__solid__density}
        "density", parameters, 1);
    DBUG("Use \'%s\' as solid density.", solid_density.name.c_str());

    // Specific heat capacity of solid.
    auto& solid_heat_c = findParameter<double>(
        mat_solid_config,
        //! \ogs_file_param_special{process__HEAT_TRANSPORT__material_property__solid__heat_capacity}
        "heat_capacity", parameters, 1);
    DBUG("Use \'%s\' as specific heat capacity of solid.",
         solid_heat_c.name.c_str());

    // Thermal conductivity of solid.
    auto& solid_thermal_c = findParameter<double>(
        mat_solid_config,
        //! \ogs_file_param_special{process__HEAT_TRANSPORT__material_property__solid__thermal_conductivity}
        "thermal_conductivity", parameters, mesh.getDimension());
    DBUG("Use \'%s\' as thermal conductivity of solid.",
         solid_thermal_c.name.c_str());

    auto const& mat_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    if (mat_ids)
    {
        INFO("The heat transport is in heterogeneous porous media.");
        return std::unique_ptr<Process>{new HeatTransportProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(secondary_variables),
            std::move(named_function_caller), solid_density, solid_heat_c,
            solid_thermal_c, *mat_ids, mat_config}};
    }
    else
    {
        INFO("The heat transport is in homogeneous porous media.");

        MeshLib::Properties dummy_property;
        auto const& dummy_property_vector =
            dummy_property.createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        return std::unique_ptr<Process>{new HeatTransportProcess{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(secondary_variables),
            std::move(named_function_caller), solid_density, solid_heat_c,
            solid_thermal_c, *dummy_property_vector, mat_config}};
    }
}

}  // end of namespace
}  // end of namespace
