/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   HeatTransportMaterialProperties.h
 *
 */

#ifndef OGS_HEAT_TRANSPORT_MATERIAL_PROPERTIES_H
#define OGS_HEAT_TRANSPORT_MATERIAL_PROPERTIES_H

#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "MaterialLib/Fluid/FluidProperty.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity;
}
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

class SpatialPosition;

namespace HeatTransport
{
class HeatTransportMaterialProperties
{
public:
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    HeatTransportMaterialProperties(
        Parameter<double> const& solid_density,
        Parameter<double> const& solid_heat_capacity,
        Parameter<double> const& thermal_conductivity,
        MeshLib::PropertyVector<int> const& material_ids,
        BaseLib::ConfigTree const& config);

    void setMaterialID(const SpatialPosition& pos);

    /**
     * \brief Compute the coefficient of the mass term
     * \param t                  Time.
     * \param pos                Position of element.
     * \param T                  Temperature value.
     * \param porosity_variable  The first variable for porosity model.
     * \param pressure           Pressure values of different phases
     * \param saturation         Saturation values of different phases
     */
    double getMassCoefficient(const double t, const SpatialPosition& pos,
                              const double T, const double porosity_variable,
                              const std::vector<double>& pressures,
                              const std::vector<double>& saturations);

    Eigen::MatrixXd getConductivity(const double t, const SpatialPosition& pos,
                                    const double p, const double T,
                                    const std::vector<double>& pressures,
                                    const std::vector<double>& saturations,
                                    const int dim) const;

private:
    struct FluidProperties
    {
        using FluidPropertyModel =
            std::unique_ptr<MaterialLib::Fluid::FluidProperty>;

        FluidProperties() = default;

        FluidProperties(FluidPropertyModel&& rho,
                        FluidPropertyModel&& c_p,
                        FluidPropertyModel&& kappa)
            : density(std::move(rho)),
              heat_capacity(std::move(c_p)),
              thermal_conductivity(std::move(kappa))
        {
        }

        ~FluidProperties() = default;

        FluidPropertyModel density;
        FluidPropertyModel heat_capacity;
        FluidPropertyModel thermal_conductivity;
    };

    /// Fluid properties for multiphase multi-component fluids
    std::vector<std::vector<std::unique_ptr<FluidProperties>>>
        _fluid_properties;
    std::vector<double> _caculated_densities;

    int _current_material_id = 0;
    /** Use porous medium models for different material zones.
     *  Material IDs must be given as mesh element properties.
     */
    MeshLib::PropertyVector<int> const& _material_ids;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    double _current_porosity = 0.0;

    // Limitation of using class Parameter for thermal conductivity:
    // The dimension of the thermal conductivity is fixed, which is not suitable
    // for problem with a hybrid of different element dimensions.
    Parameter<double> const& _thermal_conductivity;
    Parameter<double> const& _solid_density;
    Parameter<double> const& _solid_heat_capacity;
};

}  // end of namespace
}  // end of namespace
#endif /* OGS_HEAT_TRANSPORT_MATERIAL_PROPERTIES_H */
