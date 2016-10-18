/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   HeatTransportMaterialProperties.cpp
 *
 */

#include "HeatTransportMaterialProperties.h"

#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/reorderVector.h"

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

#include "MaterialLib/Fluid/Density/createFluidDensityModel.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Porosity/createPorosityModel.h"

#include "MaterialLib/Fluid/SpecificHeatCapacity/createSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/createFluidThermalConductivityModel.h"

namespace ProcessLib
{
namespace HeatTransport
{
HeatTransportMaterialProperties::HeatTransportMaterialProperties(
    Parameter<double> const& solid_density,
    Parameter<double> const& solid_heat_capacity,
    Parameter<double> const& thermal_conductivity,
    MeshLib::PropertyVector<int> const& material_ids,
    BaseLib::ConfigTree const& config)
    : _thermal_conductivity(thermal_conductivity),
      _solid_density(solid_density),
      _solid_heat_capacity(solid_heat_capacity),
      _material_ids(material_ids)
{
    DBUG("Reading material properties of heat transport process.");

    // Get fluid properties
    //! \ogs_file_param{prj__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    std::vector<int> phase_id;
    std::size_t num_fluid_comps = 0;
    //! \ogs_file_param{prj__material_property__fluid__phase}
    for (auto const& phase_fluid_conf :
         fluid_config.getConfigSubtreeList("phase"))
    {
        //! \ogs_file_attr{prj__material_property__fluid__phase__id}
        const auto id = phase_fluid_conf.getConfigAttributeOptional<int>("id");
        phase_id.push_back(*id);

        std::vector<std::unique_ptr<FluidProperties>> comp_fluid_properties;
        std::vector<int> component_id;
        //! \ogs_file_param{prj__material_property__fluid__phase__component}
        for (auto const& phase_comp_fluid_conf :
             phase_fluid_conf.getConfigSubtreeList("component"))
        {
            //! \ogs_file_attr{prj__material_property__fluid__phase__component__id}
            const auto id =
                phase_comp_fluid_conf.getConfigAttributeOptional<int>("id");
            component_id.push_back(*id);

            // FluidProperties fluid_prop;
            //! \ogs_file_param{prj__material_property__fluid__phase__component__density}
            auto const& rho_conf =
                phase_comp_fluid_conf.getConfigSubtree("density");
            auto fluid_density =
                MaterialLib::Fluid::createFluidDensityModel(rho_conf);

            //! \ogs_file_param{prj__material_property__fluid__phase__component__heat_capacity}
            auto const& c_p_conf =
                phase_comp_fluid_conf.getConfigSubtree("heat_capacity");
            auto fluid_heat_capacity =
                MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(
                    c_p_conf);

            //! \ogs_file_param{prj__material_property__fluid__phase__component__thermal_conductivity}
            auto const& beta_conf =
                phase_comp_fluid_conf.getConfigSubtree("thermal_conductivity");
            auto fluid_thermal_conductivity =
                MaterialLib::Fluid::createFluidThermalConductivityModel(
                    beta_conf);

            std::unique_ptr<FluidProperties> fluid_properties =
                std::unique_ptr<FluidProperties>(
                    new FluidProperties(std::move(fluid_density),
                                        std::move(fluid_heat_capacity),
                                        std::move(fluid_thermal_conductivity)));
            comp_fluid_properties.push_back(std::move(fluid_properties));
        }

        num_fluid_comps += comp_fluid_properties.size();

        BaseLib::reorderVector(comp_fluid_properties, component_id);
        _fluid_properties.push_back(std::move(comp_fluid_properties));
    }
    BaseLib::reorderVector(_fluid_properties, phase_id);

    _caculated_densities.resize(num_fluid_comps);

    // Get porous properties
    std::vector<int> mat_ids;
    //! \ogs_file_param{prj__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_attr{prj__material_property__porous_medium__porous_medium__id}
        auto const id = conf.getConfigAttributeOptional<int>("id");
        mat_ids.push_back(*id);
        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__porosity}
        auto const& poro_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(poro_conf);
        _porosity_models.emplace_back(std::move(n));
    }
}

void HeatTransportMaterialProperties::setMaterialID(const SpatialPosition& pos)
{
    if (_material_ids.empty())
    {
        assert(pos.getElementID().get() < _material_ids.size());
        _current_material_id = _material_ids[pos.getElementID().get()];
    }
}

double HeatTransportMaterialProperties::getMassCoefficient(
    const double t, const SpatialPosition& pos, const double T,
    const double porosity_variable, const std::vector<double>& pressures,
    const std::vector<double>& saturations)
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;

    double cp_rho_fluid = 0.;
    std::size_t phase_id = 0;
    std::size_t comp_id = 0;
    for (auto const& fluid_phase_prop : _fluid_properties)
    {
        double cp_rho_phase = 0.0;
        const double p = pressures[phase_id];
        for (auto const& fluid_comp_prop : fluid_phase_prop)
        {
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::pl)] = p;
            const double rho_c = fluid_comp_prop->density->getValue(vars);

            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::rho)] = rho_c;
            const double cp_c = fluid_comp_prop->heat_capacity->getValue(vars);
            cp_rho_phase += rho_c * cp_c;

            _caculated_densities[comp_id] = rho_c;
            comp_id++;
        }

        cp_rho_fluid += saturations[phase_id] * cp_rho_phase;
        phase_id++;
    }

    const double rho_s = _solid_density(t, pos)[0];
    const double cp_s = _solid_heat_capacity(t, pos)[0];

    // Save the porosity in the member data to be used in getConductivity
    const double _current_porosity =
        _porosity_models[_current_material_id]->getValue(porosity_variable, T);

    return _current_porosity * cp_rho_fluid +
           (1. - _current_porosity) * rho_s * cp_s;
}

Eigen::MatrixXd HeatTransportMaterialProperties::getConductivity(
    const double t, const SpatialPosition& pos, const double p, const double T,
    const std::vector<double>& pressures,
    const std::vector<double>& saturations, const int dim) const
{
    // Conductivity of liquid
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;

    double thermal_conductivity_fluid = 0.;
    std::size_t phase_id = 0;
    std::size_t comp_id = 0;
    for (auto const& fluid_phase_prop : _fluid_properties)
    {
        double k_phase = 0.0;
        const double p = pressures[phase_id];
        for (auto const& fluid_comp_prop : fluid_phase_prop)
        {
            const double rho_c = _caculated_densities[comp_id];

            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::rho)] = rho_c;
            k_phase += fluid_comp_prop->thermal_conductivity->getValue(vars);
            comp_id++;
        }

        thermal_conductivity_fluid += saturations[phase_id] * k_phase;
        phase_id++;
    }

    // Conductivity of solid
    auto k_elems_solid = _thermal_conductivity(t, pos);
    assert(k_elems_solid.size() == dim * dim);

    // Conductivity of porous medium
    for (std::size_t i = 0; i < dim; i++)
        k_elems_solid[i] =
            _current_porosity * thermal_conductivity_fluid +
            (1. - _current_porosity) * k_elems_solid[i * dim + i];

    return MathLib::toMatrix(k_elems_solid, dim, dim);
}

}  // end of namespace
}  // end of namespace
