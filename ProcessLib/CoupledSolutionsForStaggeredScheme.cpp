/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CoupledSolutionsForStaggeredScheme.cpp
 *
 * Created on November 7, 2016,
 *
 */

#include "CoupledSolutionsForStaggeredScheme.h"

#include "MathLib/LinAlg/LinAlg.h"
#include "Process.h"

namespace ProcessLib
{
CoupledSolutionsForStaggeredScheme::CoupledSolutionsForStaggeredScheme(
    std::reference_wrapper<GlobalVector const> const& coupled_xs_,
    const double dt_)
    : coupled_xs(coupled_xs_), dt(dt_)
{
    for (auto const& coupled_x : coupled_xs)
    {
        MathLib::LinAlg::setLocalAccessibleVector(coupled_x.get());
    }
}

}  // end of ProcessLib
