/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "materiallib_export.h"

#include "BaseLib/ConfigTree.h"
#include "Reaction.h"
#include "Adsorption.h"

namespace ProcessLib
{
template<typename>
class TESFEMReactionAdaptorCaOH2;
}

namespace Adsorption
{

class ReactionCaOH2 final : public Reaction
{
public:
    explicit ReactionCaOH2(BaseLib::ConfigTree const& conf) :
        //! \ogs_file_param{material__adsorption__reaction__CaOH2__ode_solver_config}
        _ode_solver_config{conf.getConfigSubtree("ode_solver_config")}
    {}

    double getEnthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override;

    double getReactionRate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override;

    const BaseLib::ConfigTree& getOdeSolverConfig() const { return _ode_solver_config; }

    // TODO merge with getReactionRate() above
    double getReactionRate(double const solid_density);

    void updateParam(double T_solid, double p_gas, double x_react,
                     double rho_s_initial);

private:
    void calculateQR();
    void setChemicalEquilibrium();
    double CaHydration();

    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    double _rho_s = nan;    //!< solid phase density
    double _p_gas = nan;    //!< gas phase pressure in unit bar
    double _p_r_g = nan;    //!< pressure of H2O on gas phase
    double _p_eq = 1.0;     //!< equilibrium pressure in bar
    double _T_eq = nan;     //!< equilibrium temperature
    double _T_s = nan;      //!< solid phase temperature
    double _qR = nan;       //!< rate of solid density change
    double _x_react = nan;  //!< mass fraction of water in gas phase
    double _X_D =
        nan;  //!< mass fraction of dehydration (CaO) in the solid phase
    double _X_H = nan;  //!< mass fraction of hydration in the solid phase

    const BaseLib::ConfigTree _ode_solver_config;

    template<typename>
    friend class ProcessLib::TESFEMReactionAdaptorCaOH2;

public:
    static MATERIALLIB_EXPORT constexpr double rho_low =
        1656.0;  //!< lower density limit
    static MATERIALLIB_EXPORT constexpr double rho_up =
        2200.0;  //!< upper density limit
};

}  // namespace Adsorption
