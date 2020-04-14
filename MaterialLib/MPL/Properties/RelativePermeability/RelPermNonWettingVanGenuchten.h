/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#pragma once

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;

/**
 *  van Genuchten relative permeability function for non-wetting phase
 *
 *   \f[k_{rel}^n= (1 - S_e)^{1/3} (1 - S_e^{1/m})^{2m}\f]
 *   with
 *   \f[S_e=\frac{S^w-S_r}{S^w_{\mbox{max}}-S^w_r}\f]
 *   where
 *    \f{eqnarray*}{
 *       &S^w_r&            \mbox{residual saturation of wetting phase,}\\
 *       &S^w_{\mbox{max}}& \mbox{maximum saturation of wetting phase,}\\
 *       &m\, \in (0, 1) &    \mbox{ exponent.}\\
 *    \f}
 */
class RelPermNonWettingVanGenuchten final : public Property
{
public:
    /**
    * @param Snr       Residual saturation of the non-wetting phase,
    *                  \f$ S^n_r \f$
    * @param Snmax     Maximum saturation  of the non-wetting phase,
    *                  \f$ S^n_{\mbox{max}} \f$
    * @param m         Exponent, \f$ m \in [0,1]\f$
    * @param krel_min  Minimum relative permeability,
    *                  \f$ k_{rel}^{\mbox{min}}\f$
    */
    RelPermNonWettingVanGenuchten(const double Snr, const double Snmax,
                                  const double m, const double krel_min);
    /// This method assigns a pointer to the material object that is the owner
    /// of this property
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (std::holds_alternative<Phase*>(scale_pointer))
        {
            _phase = std::get<Phase*>(scale_pointer);
        }
        else
        {
            OGS_FATAL(
                "The property 'RelPermNonWettingVanGenuchten' is "
                "implemented on the 'phase' scale only.");
        }
    }

    /// \return \f$k_{rel}^n \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// \return \f$ \dfrac{\partial k_{rel}^n}{\partial S^w} \f$.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    Phase* _phase = nullptr;

    /// A small number for an offset to set the bound of S, the saturation, such
    /// that S in  [Sr+_minor_offset, Smax-_minor_offset].
    const double _minor_offset = std::numeric_limits<double>::epsilon();

    const double _saturation_r;    ///< Residual saturation of wetting phase.
    const double _saturation_max;  ///< Maximum saturation of wetting phase.
    const double _m;               ///< Exponent \f$ m \f$.
    const double _krel_min;        ///< Minimum relative permeability
};
}  // namespace MaterialPropertyLib
