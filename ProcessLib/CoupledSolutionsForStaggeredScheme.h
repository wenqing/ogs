/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CouplingSolutionsForStaggeredScheme.h
 *
 * Created on November 7, 2016, 12:14 PM
 */

#pragma once

#include <functional>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
class Process;

/**
 *  A struct to keep the references of the current solutions of the equations of
 *  the coupled processes.
 *
 *  During staggered coupling iteration, an instance of this struct is created
 *  and passed through interfaces to global and local assemblers for each
 *  process.
 */
struct CoupledSolutionsForStaggeredScheme
{
    CoupledSolutionsForStaggeredScheme(
        std::vector<std::reference_wrapper<GlobalVector const>> const&
            coupled_xs_,
        const double dt_);

    /// References to the current solutions of the coupled processes.
    std::vector<std::reference_wrapper<GlobalVector const>> const& coupled_xs;

    const double dt;  ///< Time step size.
};

}  // end of ProcessLib
