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

#include "createFluidThermalConductivityModel.h"

#include "BaseLib/Error.h"
#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"

#include "MaterialLib/Fluid/ThermalConductivity/ThermalConductivityWaterIAPWS.h"

namespace MaterialLib
{
namespace Fluid
{
std::unique_ptr<FluidProperty> createFluidThermalConductivityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__thermal_conductivity__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
        return std::unique_ptr<FluidProperty>(new ConstantFluidProperty(
            //! \ogs_file_param{material__fluid__thermal_conductivity__Constant__value}
            config.getConfigParameter<double>("value")));
    else if (type == "WaterIAPWS")
        return std::unique_ptr<FluidProperty>(
            new ThermalConductivityWaterIAPWS());
    // TODO: add more models
    else
    {
        OGS_FATAL(
            "The viscosity type %s is unavailable.\n"
            "The available type is \n\tConstant \n\tWaterIAPWS",
            type.data());
    }
}

}  // end namespace
}  // end namespace
