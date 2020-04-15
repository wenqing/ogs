/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 27, 2020, 4:01 PM
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <random>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/CreateRelPermNonWettingVanGenuchten.h"
#include "MaterialLib/MPL/Properties/RelativePermeability/RelPermNonWettingVanGenuchten.h"

#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

std::unique_ptr<MPL::Property> createTestProperty(
    const char xml[],
    std::function<std::unique_ptr<MaterialPropertyLib::Property>(
        BaseLib::ConfigTree const& config)>
        createProperty)
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("property");
    // Parsing the property name:
    auto const property_name =
        sub_config.getConfigParameter<std::string>("name");

    return createProperty(sub_config);
}

TEST(TestNonWettingPhaseRelativePermeability, RelPermNonWettingVanGenuchten)
{
    const char xml[] =
        "<property>"
        "   <name>relative_permeability</name>"
        "   <type>RelativePermeabilityNonWettingVanGenuchten</type>"
        "   <residual_saturation>0.1</residual_saturation>"
        "   <maximum_saturation>1.0</maximum_saturation>"
        "   <exponent>0.5</exponent>"
        "   "
        "<minimum_relative_permeability>1.e-9</minimum_relative_permeability>"
        "</property>";

    std::unique_ptr<MPL::Property> const perm_ptr =
        createTestProperty(xml, MPL::createRelPermNonWettingVanGenuchten);
    MPL::Property const& perm_model = *perm_ptr;

    std::vector<double> const S = {0.2, 0.33, 0.45, 0.52, 0.6, 0.85};
    std::vector<double> const expected_krel = {
        0.87422700239237161, 0.74331414436457388, 0.59527539448807487,
        0.49976666464188485, 0.38520070797257489, 0.041219134248319585};

    const double perturbation = 1.e-9;
    for (std::size_t i = 0; i < S.size(); i++)
    {
        MPL::VariableArray variable_array;
        ParameterLib::SpatialPosition const pos;
        double const t = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i];

        const double k_rel_i =
            perm_model.template value<double>(variable_array, pos, t, dt);
        ASSERT_NEAR(expected_krel[i], k_rel_i, 1.e-9);

        const double dkrel_dS = perm_model.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::liquid_saturation,
            pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i] + perturbation;
        const double k_rel_i_1 =
            perm_model.template value<double>(variable_array, pos, t, dt);

        // Compare the derivative with numerical one.
        ASSERT_NEAR(dkrel_dS, (k_rel_i_1 - k_rel_i) / perturbation, 1.e-6);
    }
}
