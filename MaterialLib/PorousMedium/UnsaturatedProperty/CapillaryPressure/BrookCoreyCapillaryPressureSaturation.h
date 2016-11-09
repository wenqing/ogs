/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   BrookCoreyCapillaryPressureSaturation.h
 *
 *  Created on November 1, 2016, 9:45 AM
 */

#ifndef OGS_BROOK_COREY_CAPILLARY_PRESSURE_SATURATION_H
#define OGS_BROOK_COREY_CAPILLARY_PRESSURE_SATURATION_H

#include <limits>  // std::numeric_limits

#include "CapillaryPressureSaturation.h"

namespace MaterialLib
{
namespace PorousMedium
{
/**
 *   \brief Brook-Corey capillary pressure saturation model
 *
 *   \f[p_c=p_b S_e^{-1/m}\f]
 *   with
 *   \f[S_e=\dfrac{S-S_r}{S_{\mbox{max}}-S_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &p_b&            \mbox{ entry pressure,}\\
 *       &S_r&            \mbox{ residual saturation,}\\
 *       &S_{\mbox{max}}& \mbox{ maximum saturation,}\\
 *       &m(>=1) &        \mbox{ exponent.}\\
 *    \f}
 *
 *    Note:
 *     \f[m=1/(1-n)\f]
 */
class BrookCoreyCapillaryPressureSaturation final
    : public CapillaryPressureSaturation
{
public:
    /**
     * @param pb     Entry pressure, \f$ p_b \f$
     * @param Sr     Residual saturation, \f$ S_r \f$
     * @param Smax   Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param m      Exponent, \f$ m \f$
     * @param Pc_max Maximum capillary pressure, \f$ P_c^{\mbox{max}}\f$
     */
    BrookCoreyCapillaryPressureSaturation(const double pb, const double Sr,
                                          const double Smax, const double m,
                                          const double Pc_max)
        : _pb(pb), _Sr(Sr), _Smax(Smax), _mm(m), _Pc_max(Pc_max)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Brook-Corey capillary pressure saturation model.";
    }

    /// Get capillary pressure.
    double getCapillaryPressure(const double saturation) const override;

    /// Get capillary pressure.
    double getSaturation(const double capillary_pressure) const override;

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override;

private:
    const double _pb;      ///< Entry pressure.
    const double _Sr;      ///< Residual saturation.
    const double _Smax;    ///< Maximum saturation.
    const double _mm;      ///< Exponent (<=1.0), n=1/(1-mm).
    const double _Pc_max;  ///< Maximum capillaray pressure

    /** A small number for an offset:
     *  1. to set the bound of S, the saturation, such that
     *     S in  [_Sr+_minor_offset, _Smax-_minor_offset]
     *  2. to set the bound of Pc, the capillary pressure, such that
     *     Pc in [_minor_offset, _Pc_max]
     */
    const double _minor_offset = std::numeric_limits<double>::epsilon();
};

}  // end namespace
}  // end namespace
#endif /* OGS_BROOK_COREY_CAPILLARY_PRESSURE_SATURATION_H */
