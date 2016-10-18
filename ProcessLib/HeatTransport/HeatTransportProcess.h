/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   HeatTransportProcess.h
 *
 */

#ifndef OGS_HEAT_TRANSPORT_PROCESS_H
#define OGS_HEAT_TRANSPORT_PROCESS_H

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"

#include "HeatTransportMaterialProperties.h"
#include "HeatTransportLocalAssembler.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace HeatTransport
{
/**
 * \brief A class to simulate the heat transport process in porous media
 * described
 * by
 *
 * \f[
 *     (n \rho^l c_p^l + n \rho^s c_p^s)  \frac{\partial T}{\partial t}
 *      -\nabla (\mathbf K \nabla T) + \mathbf V \cdot \nabla T= Q
 * \f]
 * where
 *    \f{eqnarray*}{
 *       &T: &         \mbox{Temperature,}\\
 *       &\rho^l:&     \mbox{density of liquid,}\\
 *       &c_p^l:&      \mbox{specific capacity of liquid,}\\
 *       &\rho^s:&     \mbox{density of solid,}\\
 *       &c_p^s:&      \mbox{specific capacity of solid,}\\
 *       &n:&          \mbox{porosity,}\\
 *       &\mathbf K:&  \mbox{thermal conductivity,}\\
 *       &\mathbf V:&  \mbox{Fluid velocity}\\
 *    \f}
 */
class HeatTransportProcess final : public Process
{
public:
    HeatTransportProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        Parameter<double> const& solid_density,
        Parameter<double> const& solid_heat_capacity,
        Parameter<double> const& thermal_conductivity,
        MeshLib::PropertyVector<int> const& material_ids,
        BaseLib::ConfigTree const& config);

    bool isLinear() const override { return true; }
private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    HeatTransportMaterialProperties _material_properties;

    std::vector<std::unique_ptr<HeatTransportLocalAssemblerInterface>>
        _local_assemblers;
};

}  // end of namespace
}  // end of namespace

#endif  // OGS_HEAT_TRANSPORT_PROCESS_H
