/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   HTLocalAssemblerInterface.h
 *  Created on October 11, 2017, 1:35 PM
 */

#pragma once

#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
namespace HT
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType N_,
                         GlobalDimNodalMatrixType dNdx_,
                         double const& integration_weight_)
        : N(std::move(N_)),
          dNdx(std::move(dNdx_)),
          integration_weight(integration_weight_)
    {
    }

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class HTLocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                  public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;
};

}  // namespace HT
}  // namespace ProcessLib