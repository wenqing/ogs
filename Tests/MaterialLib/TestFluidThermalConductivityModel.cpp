/**
 *  \brief Heat conductivity model of water according to
 *         the IAPWS Industrial Formulation 1997
 *         http://www.iapws.org/relguide/IF97-Rev.html
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file TestFluidThermalConductivityModel.cpp
 *
 */

#include <gtest/gtest.h>

#include <memory>

#include "TestTools.h"

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "MaterialLib/Fluid/ThermalConductivity/createFluidThermalConductivityModel.h"
#include "MaterialLib/Fluid/Density/LiquidDensity.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperty> createFluidThermalConductivityModel(
    const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("thermal_conductivity");
    return MaterialLib::Fluid::createFluidThermalConductivityModel(sub_config);
}

TEST(Material, checkConstantFluidThermalConductivity)
{
    const char xml[] =
        "<thermal_conductivity>"
        "   <type>Constant</type>"
        "   <value> .45 </value> "
        "</thermal_conductivity>";
    const auto lambda = createFluidThermalConductivityModel(xml);

    ArrayType dummy;
    ASSERT_EQ(.45, lambda->getValue(dummy));
}

TEST(Material, checkFluidThermalConductivityWaterIAPWS)
{
    std::array<double, 5> const parameters{
        {2.0e-4, 999.8, 273.15, 1.e+5, 2.15e9}};
    MaterialLib::Fluid::LiquidDensity density(parameters);

    ArrayType vars;
    vars[static_cast<int>(PropertyVariableType::T)] =
        MaterialLib::PhysicalConstant::CelsiusZeroInKelvin + 20.0;
    vars[static_cast<int>(PropertyVariableType::pl)] = 1.e+5;
    const double current_rho = density.getValue(vars);
    vars[static_cast<int>(PropertyVariableType::rho)] = current_rho;

    const char xml[] =
        "<thermal_conductivity>"
        "   <type>WaterIAPWS</type>"
        "</thermal_conductivity>";
    const auto beta = createFluidThermalConductivityModel(xml);

    // Thermal conductivity of water at T=20 °C and P = 1 bar
    // should be close to 0.6.
    ASSERT_NEAR(0.598, beta->getValue(vars), 1.e-3);
}
