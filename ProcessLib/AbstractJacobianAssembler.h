/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "BaseLib/Error.h"

namespace ProcessLib
{
class LocalAssemblerInterface;

//! Base class for Jacobian assemblers.
class AbstractJacobianAssembler
{
public:
    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$.
    virtual void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                                      double const t, double const dt,
                                      std::vector<double> const& local_x,
                                      std::vector<double> const& local_x_prev,
                                      std::vector<double>& local_b_data,
                                      std::vector<double>& local_Jac_data) = 0;

    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$ with coupling.
    virtual void assembleWithJacobianForStaggeredScheme(
        LocalAssemblerInterface& /*local_assembler*/, double const /*t*/,
        double const /*dt*/, Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/, int const /*process_id*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
    {
        // TODO make pure virtual.
        OGS_FATAL("not implemented.");
    }

    virtual std::unique_ptr<AbstractJacobianAssembler> copy() const = 0;

    virtual ~AbstractJacobianAssembler() = default;
};

}  // namespace ProcessLib
