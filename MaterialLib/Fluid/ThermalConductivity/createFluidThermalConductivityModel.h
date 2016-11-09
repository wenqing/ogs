/**
 *  \brief A function for creating a thermal conductivity model for fluid
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file createFluidThermalConductivityModel.h
 *
 */

#ifndef OGS_CREATE_FLUID_THERMAL_CONDUCTIVITY_MODEL_H
#define OGS_CREATE_FLUID_THERMAL_CONDUCTIVITY_MODEL_H

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace Fluid
{
class FluidProperty;

/**
 *  Create a viscosity model
 *  \param config  ConfigTree object has a tag of <fluid_thermal_conductivity>
 */
std::unique_ptr<FluidProperty> createFluidThermalConductivityModel(
    BaseLib::ConfigTree const& config);

}  // end namespace
}  // end namespace

#endif  // OGS_CREATE_FLUID_THERMAL_CONDUCTIVITY_MODEL_H
