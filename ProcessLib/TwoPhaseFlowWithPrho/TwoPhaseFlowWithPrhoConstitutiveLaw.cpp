/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TwoPhaseFlowWithPrhoConstitutiveLaw.h"

#include <utility>

#include "NumLib/NewtonRaphson.h"

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/MPL/Property.h"

using MaterialLib::PhysicalConstant::MolarMass::H2;
using MaterialLib::PhysicalConstant::IdealGasConstant;
using MaterialLib::PhysicalConstant::HenryConstant::HenryConstantH2;
namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
bool TwoPhaseFlowWithPrhoConstitutiveLaw::compute(
    double const t,
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
    double& dxm_dX)
{
    {  // Local Newton solver
        using LocalJacobianMatrix =
            Eigen::Matrix<double, 2, 2, Eigen::RowMajor>;
        using LocalResidualVector = Eigen::Matrix<double, 2, 1>;
        using LocalUnknownVector = Eigen::Matrix<double, 2, 1>;
        LocalJacobianMatrix J_loc;

        Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(2);
        auto const update_residual = [&](LocalResidualVector& residual) {
            calculateResidual(t, dt, pos, capillary_pressure, variables, pl, X,
                              T, Sw, X_m, residual);
        };

        auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
            calculateJacobian(t, dt, pos, capillary_pressure, variables, pl, X,
                              T, jacobian, Sw,
                              X_m);  // for solution dependent Jacobians
        };

        auto const update_solution = [&](LocalUnknownVector const& increment) {
            // increment solution vectors
            Sw += increment[0];
            X_m += increment[1];
        };

        // TODO Make the following choice of maximum iterations and convergence
        // criteria available from the input file configuration. See Ehlers
        // material model implementation for the example.
        const int maximum_iterations(20);
        const double tolerance(1.e-14);

        auto newton_solver = NumLib::NewtonRaphson<
            decltype(linear_solver), LocalJacobianMatrix,
            decltype(update_jacobian), LocalResidualVector,
            decltype(update_residual), decltype(update_solution)>(
            linear_solver, update_jacobian, update_residual, update_solution,
            {maximum_iterations, tolerance});

        auto const success_iterations = newton_solver.solve(J_loc);

        if (!success_iterations)
        {
            return false;
        }
    }

    variables[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

    const double pg =
        pl + capillary_pressure.value<double>(variables, pos, t, dt);
    double const dPC_dSw = capillary_pressure.dValue<double>(
        variables, MaterialPropertyLib::Variable::liquid_saturation, pos, t,
        dt);

    dsw_dpg = calculatedSwdP(pg, Sw, dPC_dSw, X_m, T);
    dsw_dX = calculatedSwdX(pg, Sw, dPC_dSw, X_m, T);
    dxm_dpg = calculatedXmdP(pg, Sw, dPC_dSw, X_m, dsw_dpg);
    dxm_dX = calculatedXmdX(pg, Sw, dPC_dSw, X_m, dsw_dX);
    return true;
}
void TwoPhaseFlowWithPrhoConstitutiveLaw::calculateResidual(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    MaterialPropertyLib::Property const& capillary_pressure,
    MaterialPropertyLib::VariableArray& variables, double const pl,
    double const X, double const T, double Sw, double rho_h2_wet,
    ResidualVector& res)
{
    variables[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

    const double pg =
        pl + capillary_pressure.value<double>(variables, pos, t, dt);
    const double rho_h2_nonwet = pg * H2 / IdealGasConstant / T;

    // calculating residual
    res(0) = calculateEquilibiumRhoWetLight(pg, Sw, rho_h2_wet);
    res(1) = calculateSaturation(X, Sw, rho_h2_wet, rho_h2_nonwet);
}

void TwoPhaseFlowWithPrhoConstitutiveLaw::calculateJacobian(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    MaterialPropertyLib::Property const& capillary_pressure,
    MaterialPropertyLib::VariableArray& variables, double const pl,
    double const /*X*/, double const T, JacobianMatrix& Jac, double Sw,
    double rho_h2_wet)
{
    variables[static_cast<int>(
        MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

    const double pg =
        pl + capillary_pressure.value<double>(variables, pos, t, dt);
    const double rho_h2_nonwet = pg * H2 / IdealGasConstant / T;
    double const rho_equili_h2_wet = pg * HenryConstantH2 * H2;
    double const dPC_dSw = capillary_pressure.dValue<double>(
        variables, MaterialPropertyLib::Variable::liquid_saturation, pos, t,
        dt);
    double const drhoh2wet_dpg = HenryConstantH2 * H2;
    Jac.setZero();
    if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
    {
        Jac(0, 0) = -1;
    }
    else
    {
        Jac(0, 0) = drhoh2wet_dpg * dPC_dSw;
        Jac(0, 1) = -1;
    }

    Jac(1, 0) = rho_h2_nonwet - rho_h2_wet;
    Jac(1, 1) = -Sw;
}

double TwoPhaseFlowWithPrhoConstitutiveLaw::calculateEquilibiumRhoWetLight(
    double const pg, double const Sw, double const rho_wet_h2) const
{
    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
    return std::min(1 - Sw, rho_equilibrium_wet_h2 - rho_wet_h2);
}

double TwoPhaseFlowWithPrhoConstitutiveLaw::calculateSaturation(
    double const X, double const Sw, double const rho_wet_h2,
    double const rho_nonwet_h2) const
{
    return X - (Sw * rho_wet_h2 + (1 - Sw) * rho_nonwet_h2);
}

double TwoPhaseFlowWithPrhoConstitutiveLaw::calculatedSwdP(double const pg,
                                                           double const Sw,
                                                           const double dPC_dSw,
                                                           double rho_wet_h2,
                                                           double const T) const
{
    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
    if ((1 - Sw) < (rho_equilibrium_wet_h2 - rho_wet_h2))
    {
        return 0.0;
    }
    double const drhoh2wet_dpg = HenryConstantH2 * H2;
    double const drhoh2nonwet_dpg = H2 / IdealGasConstant / T;
    double const alpha =
        ((drhoh2nonwet_dpg - drhoh2wet_dpg) * (1 - Sw) + drhoh2wet_dpg);
    double const beta = (drhoh2nonwet_dpg - drhoh2wet_dpg) *
                        pg;  // NOTE here should be PG^h, but we ignore vapor
    return alpha / (beta - alpha * dPC_dSw);
}

double TwoPhaseFlowWithPrhoConstitutiveLaw::calculatedSwdX(
    double const pg, const double Sw, double const dPC_dSw,
    const double rho_wet_h2, double const T) const
{
    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
    if ((1 - Sw) < (rho_equilibrium_wet_h2 - rho_wet_h2))
    {
        return 0.0;
    }
    double const drhoh2wet_dpg = HenryConstantH2 * H2;
    double const drhoh2nonwet_dpg = H2 / IdealGasConstant / T;
    double const alpha =
        ((drhoh2nonwet_dpg - drhoh2wet_dpg) * (1 - Sw) + drhoh2wet_dpg);
    double const beta = (drhoh2nonwet_dpg - drhoh2wet_dpg) *
                        pg;  // NOTE here should be PG^h, but we ignore vapor
    return -1 / (beta - alpha * dPC_dSw);
}

double TwoPhaseFlowWithPrhoConstitutiveLaw::calculatedXmdX(double const pg,
                                                           double const Sw,
                                                           double const dPC_dSw,
                                                           double rho_wet_h2,
                                                           double dSwdX) const
{
    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
    if ((1 - Sw) < (rho_equilibrium_wet_h2 - rho_wet_h2))
    {
        return 1.0;
    }
    return HenryConstantH2 * H2 * dPC_dSw * dSwdX;
}

double TwoPhaseFlowWithPrhoConstitutiveLaw::calculatedXmdP(double const pg,
                                                           double const Sw,
                                                           double const dPC_dSw,
                                                           double rho_wet_h2,
                                                           double dSwdP) const
{
    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
    if ((1 - Sw) < (rho_equilibrium_wet_h2 - rho_wet_h2))
    {
        return 0.0;
    }
    return HenryConstantH2 * H2 * (1 + dPC_dSw * dSwdP);
}

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
