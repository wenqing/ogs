/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionResidual.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
ConvergenceCriterionResidual::ConvergenceCriterionResidual(
    std::optional<double>&& absolute_tolerance,
    std::optional<double>&& relative_tolerance,
    const MathLib::VecNormType norm_type)
    : ConvergenceCriterion(norm_type),
      _abstol(std::move(absolute_tolerance)),
      _reltol(std::move(relative_tolerance))
{
    if ((!_abstol) && (!_reltol))
    {
        OGS_FATAL(
            "At least one of absolute or relative tolerance has to be "
            "specified.");
    }
}

void ConvergenceCriterionResidual::checkDeltaX(
    const GlobalVector& minus_delta_x, GlobalVector const& x)
{
    auto error_dx = MathLib::LinAlg::norm(minus_delta_x, _norm_type);
    auto norm_x = MathLib::LinAlg::norm(x, _norm_type);

    INFO("Convergence criterion: |dx|={:.4e}, |x|={:.4e}, |dx|/|x|={:.4e}",
         error_dx, norm_x,
         (norm_x == 0. ? std::numeric_limits<double>::quiet_NaN()
                       : (error_dx / norm_x)));
}

void ConvergenceCriterionResidual::checkResidual(const GlobalVector& residual)
{
    auto norm_res = MathLib::LinAlg::norm(residual, _norm_type);

    if (_is_first_iteration)
    {
        INFO("Convergence criterion: |r0|={:.4e}", norm_res);
        _residual_norm_0 = norm_res;
    }
    else
    {
        _residual_norm_0 =
            (_residual_norm_0 < std::numeric_limits<double>::epsilon())
                ? norm_res
                : _residual_norm_0;
        if (_residual_norm_0 < std::numeric_limits<double>::epsilon())
        {
            INFO("Convergence criterion: |r|={:.4e} |r0|={:.4e}", norm_res,
                 _residual_norm_0);
        }
        else
        {
            INFO(
                "Convergence criterion: |r|={:.4e} |r0|={:.4e} |r|/|r0|={:.4e}",
                norm_res, _residual_norm_0,
                (_residual_norm_0 == 0.
                     ? std::numeric_limits<double>::quiet_NaN()
                     : (norm_res / _residual_norm_0)));
        }
    }

    bool satisfied_abs = false;
    bool satisfied_rel = false;

    if (_abstol)
    {
        satisfied_abs = norm_res < *_abstol;
    }
    if (_reltol && !_is_first_iteration)
    {
        satisfied_rel =
            checkRelativeTolerance(*_reltol, norm_res, _residual_norm_0);
    }

    _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
}

std::unique_ptr<ConvergenceCriterionResidual>
createConvergenceCriterionResidual(const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__type}
    config.checkConfigParameter("type", "Residual");

    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__Residual__abstol}
    auto abstol = config.getConfigParameterOptional<double>("abstol");
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__Residual__reltol}
    auto reltol = config.getConfigParameterOptional<double>("reltol");
    auto const norm_type_str =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__Residual__norm_type}
        config.getConfigParameter<std::string>("norm_type");
    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
    {
        OGS_FATAL("Unknown vector norm type `{:s}'.", norm_type_str);
    }

    return std::make_unique<ConvergenceCriterionResidual>(
        std::move(abstol), std::move(reltol), norm_type);
}

}  // namespace NumLib
