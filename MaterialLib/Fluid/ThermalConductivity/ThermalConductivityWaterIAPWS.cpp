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
 *   \file ThermalConductivityWaterIAPWS.cpp
 *
 */
#include "ThermalConductivityWaterIAPWS.h"

#include <cmath>

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace Fluid
{
double ThermalConductivityWaterIAPWS::getValue(const ArrayType& var_vals) const
{
    const double T =
        var_vals[static_cast<int>(PropertyVariableType::T)] / 647.096;
    const double rho =
        var_vals[static_cast<int>(PropertyVariableType::rho)] / 317.11;

    const double a[4] = {0.0102811, 0.0299621, 0.0156146, -0.00422464};

    double sum = 0.;
    for (int i = 0; i < 4; i++)
        sum += a[i] * MathLib::fastpow(T, i);

    const double b0 = -0.397070;
    const double b1 = 0.400302;
    const double b2 = 1.060000;
    const double B0 = -0.171587;
    const double B1 = 2.392190;
    const double lamda0 = std::sqrt(T) * sum;
    const double lamda1 =
        b0 + b1 * rho + b2 * std::exp(B0 * (rho + B1) * (rho + B1));

    const double d0 = 0.0701309;
    const double d1 = 0.0118520;
    const double d2 = 0.00169937;
    const double d3 = -1.0200;

    const double C0 = 0.642857;
    const double C1 = -4.11717;
    const double C2 = -6.17937;
    const double C3 = 0.00308976;
    const double C4 = 0.0822994;
    const double C5 = 10.0932;

    const double dT = std::fabs(T - 1) + C3;

    double S = 0.;
    if (T >= 1.0)
        S = 1 / dT;
    else
        S = C5 / std::pow(dT, 0.6);

    const double Q = 2 + (C4 / std::pow(dT, 0.6));

    const double lamda2 =
        (d0 / MathLib::fastpow(T, 10) + d1) * std::pow(rho, 1.8) *
            exp(C0 * (1 - std::pow(rho, 2.8))) +
        d2 * S * std::pow(rho, Q) *
            exp((Q / (1. + Q)) * (1 - std::pow(rho, (1. + Q)))) +
        d3 * exp(C1 * std::pow(T, 1.5) + C2 / MathLib::fastpow(rho, 5));

    return lamda0 + lamda1 + lamda2;  // lamda in [W/m/K]
}

}  // end namespace
}  // end namespace
