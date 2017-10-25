/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PhaseFieldStaggeredFEM.h"
#include "PhaseFieldStaggeredProcessData.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace PhaseFieldStaggered
{
class PhaseFieldStaggeredProcess final : public Process
{
public:
    PhaseFieldStaggeredProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        PhaseFieldStaggeredProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;

    void computeSecondaryVariableConcrete(double const t,
                                          GlobalVector const& x) override;

    void preTimestepConcreteProcess(GlobalVector const& x, const double t,
                                    const double delta_t) override;

    // Get the solution of the previous time step.
    GlobalVector* getPreviousTimeStepSolution() const override
    {
        return _x_previous_timestep.get();
    }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void postTimestepConcreteProcess(GlobalVector const& x) override;

    PhaseFieldStaggeredProcessData _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    /// Solution of the previous time step
    std::unique_ptr<GlobalVector> _x_previous_timestep = nullptr;
};

}  // namespace PhaseFieldStaggered
}  // namespace ProcessLib
