/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldSmallDeformationProcess.h"
#include "PhaseFieldSmallDeformationProcess-impl.h"

namespace ProcessLib
{
namespace PhaseFieldSmallDeformation
{
template class PhaseFieldSmallDeformationProcess<2>;
template class PhaseFieldSmallDeformationProcess<3>;

}  // namespace PhaseFieldSmallDeformation
}  // namespace ProcessLib
