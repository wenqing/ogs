/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Eigen>

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

namespace PhaseField
{
template <int DisplacementDim>
struct PhaseFieldProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::map<int, std::unique_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& kinetic_coefficient;
    ParameterLib::Parameter<double> const& solid_density;
    ParameterLib::Parameter<double>& history_field;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    bool propagating_crack = false;
    bool crack_pressure = false;
    bool constant_crack_vol = false;
    int split_method = 0;
    int secant_method = 0;
    double reg_param = -0.01;
    double pf_irrv = 0.05;
    int at_param = 2;

    double const unity_pressure = 1.0;
    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
    double injected_volume = 0.0;
    double crack_volume = 0.0;
    double pressure_n = 0.0, pressure_nm1 = 0.0;
    double crack_volume_n = 0.0, crack_volume_nm1 = 0.0;
    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    int nl_itr = 0;
};

}  // namespace PhaseField
}  // namespace ProcessLib
