/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicalPhaseFieldFEM.h"
#include "HydroMechanicalPhaseFieldProcess.h"
#include "HydroMechanicalPhaseFieldProcessData.h"

#include <cassert>

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HydroMechanicalPhaseField
{
template <int DisplacementDim>
HydroMechanicalPhaseFieldProcess<DisplacementDim>::
    HydroMechanicalPhaseFieldProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HydroMechanicalPhaseFieldProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        int const mechanics_related_process_id,
        int const phase_field_process_id,
        int const hydro_process_id)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              false),
      _process_data(std::move(process_data)),
      _mechanics_related_process_id(mechanics_related_process_id),
      _phase_field_process_id(phase_field_process_id),
      _hydro_process_id(hydro_process_id)
{
}

template <int DisplacementDim>
bool HydroMechanicalPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
HydroMechanicalPhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field process or hydro.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
HydroMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTable(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or hydro.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap&
HydroMechanicalPhaseFieldProcess<DisplacementDim>::getDOFTableByProcessID(
    const int process_id) const
{
    if (process_id == _mechanics_related_process_id)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield or hydro.
    return *_local_to_global_index_map_single_component;
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);

    // For displacement equation.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    getProcessVariables(_mechanics_related_process_id)[0]
                        .get()
                        .getNumberOfComponents(),
                    [&]() { return *_mesh_subset_all_nodes; });

    std::vector<int> const vec_n_components{DisplacementDim};
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    // For phase field equation or the heat conduction.
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<DisplacementDim>::
    initializeConcreteProcess(NumLib::LocalToGlobalIndexMap const& dof_table,
                              MeshLib::Mesh const& mesh,
                              unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, HydroMechanicalPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _mechanics_related_process_id, _phase_field_process_id,
        _hydro_process_id);

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &HydroMechanicalPhaseFieldLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &HydroMechanicalPhaseFieldLocalAssemblerInterface::
                             getIntPtEpsilon));

    /*    _secondary_variables.addSecondaryVariable(
            "velocity",
            makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                             _local_assemblers,
                             &HydroMechanicalPhaseFieldLocalAssemblerInterface::
                                 getIntPtDarcyVelocity));
                                 */

    _process_data.ele_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "damage", MeshLib::MeshItemType::Cell,
        1);

    _process_data.ele_u_dot_grad_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "u_dot_grad_d", MeshLib::MeshItemType::Cell,
        1);

    _process_data.width = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "width", MeshLib::MeshItemType::Cell,
        1);

    _process_data.width_prev = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "width_prev", MeshLib::MeshItemType::Cell,
        1);

    _process_data.cum_grad_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "cum_grad_d", MeshLib::MeshItemType::Cell,
        1);
    _process_data.ele_grad_d = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "grad_damage",
        MeshLib::MeshItemType::Cell, DisplacementDim);
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<
    DisplacementDim>::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for the equations of temperature-deformation.
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_mechanics_related_process_id),
        _mechanics_related_process_id);
    // for the phase field
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_phase_field_process_id),
        _phase_field_process_id);
    // for heat conduction
    initializeProcessBoundaryConditionsAndSourceTerms(
        getDOFTableByProcessID(_hydro_process_id), _hydro_process_id);
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for HydroMechanicalPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the staggered scheme
    if (_coupled_solutions->process_id == _mechanics_related_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "poroelastic-deformation in "
            "HydroMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }
    else if (_coupled_solutions->process_id == _phase_field_process_id)
    {
        DBUG(
            "Assemble the Jacobian equations of z"
            "phase field in "
            "HydroMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of "
            "fluid flow in "
            "HydroMechanicalPhaseFieldProcess for "
            "the staggered scheme.");
    }
    dof_tables.emplace_back(getDOFTableByProcessID(_hydro_process_id));
    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<
    DisplacementDim>::preTimestepConcreteProcess(GlobalVector const& x,
                                                 double const t,
                                                 double const dt,
                                                 const int process_id)
{
    DBUG("PreTimestep HydroMechanicalPhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;
    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);

    if (process_id != _mechanics_related_process_id)
    {
        return;
    }
    GlobalExecutor::executeMemberOnDereferenced(
        &HydroMechanicalPhaseFieldLocalAssemblerInterface::preTimestep,
        _local_assemblers, getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<
    DisplacementDim>::postTimestepConcreteProcess(GlobalVector const& x,
                                                  double const /*t*/,
                                                  double const /*dt*/,
                                                  int const process_id)
{
    DBUG("PostTimestep HydroMechanicalPhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &HydroMechanicalPhaseFieldLocalAssemblerInterface::postTimestep,
        _local_assemblers, getDOFTable(process_id), x);

    /*
        _process_data.poroelastic_energy = 0.0;
        _process_data.surface_energy = 0.0;
        _process_data.pressure_work = 0.0;
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

    dof_tables.emplace_back(getDOFTableByProcessID(_hydro_process_id));
    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));

        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::computeEnergy, _local_assemblers,
     dof_tables, x, _process_data.t, _process_data.poroelastic_energy,
            _process_data.surface_energy, _process_data.pressure_work,
            _use_monolithic_scheme, _coupled_solutions);

        INFO("Elastic energy: %g Surface energy: %g Pressure work: %g ",
             _process_data.elastic_energy, _process_data.surface_energy,
             _process_data.pressure_work);
             */
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<DisplacementDim>::
    postNonLinearSolverConcreteProcess(GlobalVector const& /*x*/,
                                       const double t,
                                       const int process_id)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    dof_tables.emplace_back(getDOFTableByProcessID(_hydro_process_id));
    dof_tables.emplace_back(
        getDOFTableByProcessID(_mechanics_related_process_id));
    dof_tables.emplace_back(getDOFTableByProcessID(_phase_field_process_id));

    if (process_id == _phase_field_process_id)
    {
        DBUG("Fracture width computation before hydro process.");
        GlobalExecutor::executeMemberOnDereferenced(
            &HydroMechanicalPhaseFieldLocalAssemblerInterface::
                computeFractureNormal,
            _local_assemblers, dof_tables, _coupled_solutions);
        GlobalExecutor::executeMemberOnDereferenced(
            &HydroMechanicalPhaseFieldLocalAssemblerInterface::
                computeFractureWidth,
            _local_assemblers, dof_tables, t, _coupled_solutions, _mesh);
    }
}

template <int DisplacementDim>
void HydroMechanicalPhaseFieldProcess<DisplacementDim>::updateConstraints(
    GlobalVector& lower, GlobalVector& upper)
{
    lower.setZero();
    MathLib::LinAlg::setLocalAccessibleVector(*_x_previous_timestep);
    MathLib::LinAlg::copy(*_x_previous_timestep, upper);

    GlobalIndexType x_size = _x_previous_timestep->size();

    for (GlobalIndexType i = 0; i < x_size; i++)
        if ((*_x_previous_timestep)[i] > _process_data.pf_irrv)
            upper.set(i, 1.0);
}

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
