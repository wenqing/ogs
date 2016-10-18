/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   HeatTransportProcess.cpp
 *
 */

#include "HeatTransportProcess.h"

#include <cassert>

#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "HeatTransportLocalAssembler.h"

#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

namespace ProcessLib
{
namespace HeatTransport
{
HeatTransportProcess::HeatTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    Parameter<double> const& solid_density,
    Parameter<double> const& solid_heat_capacity,
    Parameter<double> const& thermal_conductivity,
    MeshLib::PropertyVector<int> const& material_ids,
    BaseLib::ConfigTree const& config)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _material_properties(HeatTransportMaterialProperties(
          solid_density, solid_heat_capacity, thermal_conductivity,
          material_ids, config))
{
    DBUG("Create HeatTransport process.");
}

void HeatTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<HeatTransportLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _material_properties);

    _secondary_variables.addSecondaryVariable(
        "heat_flux_x", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &HeatTransportLocalAssemblerInterface::getIntPtHeatFluxX));

    if (mesh.getDimension() > 1)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_y", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &HeatTransportLocalAssemblerInterface::getIntPtHeatFluxY));
    }
    if (mesh.getDimension() > 2)
    {
        _secondary_variables.addSecondaryVariable(
            "heat_flux_z", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &HeatTransportLocalAssemblerInterface::getIntPtHeatFluxZ));
    }
}

void HeatTransportProcess::assembleConcreteProcess(const double t,
                                                   GlobalVector const& x,
                                                   GlobalMatrix& M,
                                                   GlobalMatrix& K,
                                                   GlobalVector& b)
{
    DBUG("Assemble HeatTransportProcess.");
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b);
}

void HeatTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatTransportProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac);
}

}  // end of namespace
}  // end of namespace
