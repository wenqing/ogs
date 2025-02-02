/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <cmath>

#include "Reaction.h"

namespace Adsorption
{

class AdsorptionReaction : public Reaction
{
public:
    // TODO [CL] move those three methods to water properties class
    static double getEvaporationEnthalpy(const double T_Ads);
    static double getEquilibriumVapourPressure(const double T_Ads);

    static double getMolarFraction(double xm, double M_this, double M_other);
    static double dMolarFraction(double xm, double M_this, double M_other);

    static double getLoading(const double rho_curr, const double rho_dry);

    double getEquilibriumLoading(const double p_Ads, const double T_Ads, const double M_Ads)
    const override;

    double getEnthalpy(const double p_Ads, const double T_Ads,
                       const double M_Ads) const override;
    double getReactionRate(const double p_Ads, const double T_Ads,
                           const double M_Ads,
                           const double loading) const override;

protected:
    virtual double getAdsorbateDensity(const double T_Ads) const = 0;
    virtual double getAlphaT(const double T_Ads) const = 0;
    virtual double characteristicCurve(const double A) const = 0;
    virtual double dCharacteristicCurve(const double A) const = 0;

private:
    double getEntropy(const double T_Ads, const double A) const;
};


inline double curvePolyfrac(const double* coeffs, const double x)
{
    // TODO use Horner scheme
    return (coeffs[0] + coeffs[2] * x + coeffs[4] * std::pow(x, 2) +
            coeffs[6] * std::pow(x, 3)) /
           (1.0 + coeffs[1] * x + coeffs[3] * std::pow(x, 2) +
            coeffs[5] * std::pow(x, 3));
}

inline double dCurvePolyfrac(const double* coeffs, const double x)
{
    const double x2 = x*x;
    const double x3 = x2*x;
    const double u  = coeffs[0] + coeffs[2] * x +     coeffs[4] * x2 +     coeffs[6] * x3;
    const double du =             coeffs[2]     + 2.0*coeffs[4] * x  + 3.0*coeffs[6] * x2;
    const double v  = 1.0 + coeffs[1] * x +     coeffs[3] * x2 +     coeffs[5] * x3;
    const double dv =       coeffs[1]     + 2.0*coeffs[3] * x  + 3.0*coeffs[5] * x2;

    return (du*v - u*dv) / v / v;
}

}  // namespace Adsorption
