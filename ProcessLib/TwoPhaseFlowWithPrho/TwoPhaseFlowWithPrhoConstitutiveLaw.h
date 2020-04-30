/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>

#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
class Property;
}

namespace ParameterLib
{
class SpatialPosition;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
class TwoPhaseFlowWithPrhoConstitutiveLaw
{
public:
    TwoPhaseFlowWithPrhoConstitutiveLaw() = default;

    bool compute(double const t,
                 double const dt,
                 ParameterLib::SpatialPosition const& pos,
                 MaterialPropertyLib::Property const& capillary_pressure,
                 MaterialPropertyLib::VariableArray& variables,
                 double const pl,
                 double const X,
                 double const T,
                 double& Sw,
                 double& X_m,
                 double& dsw_dpg,
                 double& dsw_dX,
                 double& dxm_dpg,
                 double& dxm_dX);

private:
    static int const jacobian_residual_size = 2;
    using ResidualVector = Eigen::Matrix<double, jacobian_residual_size, 1>;
    using JacobianMatrix =
        Eigen::Matrix<double, jacobian_residual_size, jacobian_residual_size,
                      Eigen::RowMajor>;
    using UnknownVector = Eigen::Matrix<double, jacobian_residual_size, 1>;

    /**
     * Calculates the residual vector.
     */
    void calculateResidual(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& pos,
        MaterialPropertyLib::Property const& capillary_pressure,
        MaterialPropertyLib::VariableArray& variables, double const pl,
        double const X, double const T, double Sw, double rho_h2_wet,
        ResidualVector& res);
    /**
     * Calculates the Jacobian.
     */
    void calculateJacobian(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& pos,
        MaterialPropertyLib::Property const& capillary_pressure,
        MaterialPropertyLib::VariableArray& variables, double const pl,
        double const X, double const T, JacobianMatrix& Jac, double Sw,
        double rho_h2_wet);
    /** Complementary condition 1
     * for calculating molar fraction of light component in the liquid phase
     */
    double calculateEquilibiumRhoWetLight(double const pg, double const Sw,
                                          double const rho_wet_h2) const;
    /** Complementary condition 2
     * for calculating the saturation
     */
    double calculateSaturation(double const X, double const Sw,
                               double const rho_wet_h2,
                               double const rho_nonwet_h2) const;
    /**
     * Calculate the derivatives using the analytical way
     */
    double calculatedSwdP(double const pg, double const Sw, const double dPcdSw,
                          double rho_wet_h2, double const T) const;
    /**
     * Calculate the derivatives using the analytical way
     */
    double calculatedSwdX(double const pg, const double Sw, const double dPcdSw,
                          const double rho_wet_h2, double const T) const;
    /**
     * Calculate the derivatives using the analytical way
     */
    double calculatedXmdX(double const pg, double const Sw, const double dPcdSw,
                          double rho_wet_h2, double dSwdX) const;
    /**
     * Calculate the derivatives using the analytical way
     */
    double calculatedXmdP(double const pg, double const Sw, const double dPcdSw,
                          double rho_wet_h2, double dSwdP) const;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
