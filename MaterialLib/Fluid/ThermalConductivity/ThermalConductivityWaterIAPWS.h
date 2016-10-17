/**
 *  \brief Thermal conductivity model of water according to
 *         the IAPWS Industrial Formulation 1997
 *         http://www.iapws.org/relguide/IF97-Rev.html
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file ThermalConductivityWaterIAPWS.h
 *
 */

#ifndef OGS_THERMAL_CONDUCTIVITY_WATER_IAPWS_H
#define OGS_THERMAL_CONDUCTIVITY_WATER_IAPWS_H

#include <string>
#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
class ThermalConductivityWaterIAPWS final : public FluidProperty
{
public:
    ThermalConductivityWaterIAPWS() = default;
    ThermalConductivityWaterIAPWS(const ThermalConductivityWaterIAPWS& /*orig*/)
    {
    }

    virtual ~ThermalConductivityWaterIAPWS() = default;

    /// Get density model name.
    std::string getName() const override
    {
        return "Heat conductivity model of water: "
               "the IAPWS Industrial Formulation 1997";
    }

    /**
     *  \param var_vals  Variable values  in an array.
     *                   Its first element is temperature, and its second
     *                   element is the fluid density.
     * \return           Heat conductivity in [W/m/K]
     */
    double getValue(const ArrayType& var_vals) const override;

    /// Get the partial differential of the property value
    double getdValue(const ArrayType& /* var_vals*/,
                     const PropertyVariableType /* var */) const override
    {
        return 0.;
    }
};

}  // end namespace
}  // end namespace

#endif  // OGS_THERMAL_CONDUCTIVITY_WATER_IAPWS_H
