/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include <memory>
#include <utility>

#include "MeshLib/PropertyVector.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace PhaseFieldInSitu
{
template <int DisplacementDim>
struct PhaseFieldInSituProcessData
{
    PhaseFieldInSituProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& residual_stiffness_,
        Parameter<double> const& crack_resistance_,
        Parameter<double> const& crack_length_scale_,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        bool propagating_crack_, bool crack_pressure_, double pf_irrv_,
        int at_param_)
        : material{std::move(material_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_),
          propagating_crack(propagating_crack_),
          crack_pressure(crack_pressure_),
          pf_irrv(pf_irrv_),
          at_param(at_param_)
    {
    }

    PhaseFieldInSituProcessData(PhaseFieldInSituProcessData&& other)
        : material{std::move(other.material)},
          residual_stiffness(other.residual_stiffness),
          crack_resistance(other.crack_resistance),
          crack_length_scale(other.crack_length_scale),
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          propagating_crack(other.propagating_crack),
          crack_pressure(other.crack_pressure),
          pf_irrv(other.pf_irrv),
          at_param(other.at_param),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    PhaseFieldInSituProcessData(PhaseFieldInSituProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldInSituProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldInSituProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& residual_stiffness;
    Parameter<double> const& crack_resistance;
    Parameter<double> const& crack_length_scale;
    Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
    double const unity_pressure = 1.0;
    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
    double injected_volume = 0.0;
    double crack_volume0 = 0.0;
    double crack_volume1 = 0.0;
    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    bool propagating_crack = false;
    bool crack_pressure = false;
    double pf_irrv = 0.05;
    int at_param = 2;
};

}  // namespace PhaseFieldInSitu
}  // namespace ProcessLib
