/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
struct PhaseTransitionData
{
    double drho_GR_dp_GR = 0.;
    double drho_GR_dp_cap = 0.;
    double drho_GR_dT = 0.;

    double drho_C_GR_dp_GR = 0.;
    double drho_C_GR_dp_cap = 0.;
    double drho_C_GR_dT = 0.;

    double drho_W_GR_dp_GR = 0.;
    double drho_W_GR_dp_cap = 0.;
    double drho_W_GR_dT = 0.;

    double drho_LR_dp_GR = 0.;
    double drho_LR_dp_cap = 0.;
    double drho_LR_dT = 0.;
    double drho_LR_dp_LR = 0.;

    double drho_C_LR_dp_GR = 0.;
    double drho_C_LR_dp_cap = 0.;
    double drho_C_LR_dT = 0.;
    double drho_C_LR_dp_LR = 0.;

    double drho_W_LR_dp_GR = 0.;
    double drho_W_LR_dp_cap = 0.;
    double drho_W_LR_dT = 0.;
    double drho_W_LR_dp_LR = 0.;

    // constituent mass and molar fractions
    double xnWG = 0.;
    double xmWG = 0.;

    // mass fraction derivatives
    double dxmWG_dpGR = 0.;
    double dxmWG_dpCap = 0.;
    double dxmWG_dT = 0.;

    double dxmWL_dpGR = 0.;
    double dxmWL_dpCap = 0.;
    double dxmWL_dpLR = 0.;
    double dxmWL_dT = 0.;

    double diffusion_coefficient_vapour = 0.;
    double diffusion_coefficient_solute = 0.;

    // specific enthalpies
    double hCG = 0;
    double hWG = 0;
    double hCL = 0;
    double hWL = 0;

    double dh_G_dT = 0;
    double dh_L_dT = 0;

    // specific inner energies
    double uG = 0;
    double uL = 0;

    double du_G_dT = 0;
    double du_L_dT = 0;
    double du_G_dp_GR = 0;
    double du_L_dp_GR = 0;
    double du_L_dp_cap = 0;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M