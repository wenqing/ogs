/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerInterface.h"
#include <cassert>
#include "NumLib/DOF/DOFTableUtil.h"

#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
void LocalAssemblerInterface::assemble(double const /*t*/,
                                       std::vector<double> const& /*local_x*/,
                                       std::vector<double>& /*local_M_data*/,
                                       std::vector<double>& /*local_K_data*/,
                                       std::vector<double>& /*local_b_data*/)
{
    OGS_FATAL(
        "The assemble() function is not implemented in the local assembler.");
}

void LocalAssemblerInterface::assembleForStaggeredScheme(
    double const /*t*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    LocalCoupledSolutions const& /*coupled_solutions*/)
{
    OGS_FATAL(
        "The assembleForStaggeredScheme() function is not implemented in the "
        "local assembler.");
}

void LocalAssemblerInterface::assembleWithJacobian(
    double const /*t*/, std::vector<double> const& /*local_x*/,
    std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
    const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL(
        "The assembleWithJacobian() function is not implemented in the local "
        "assembler.");
}

void LocalAssemblerInterface::assembleWithJacobianForStaggeredScheme(
    double const /*t*/, std::vector<double> const& /*local_xdot*/,
    const double /*dxdot_dx*/, const double /*dx_dx*/,
    std::vector<double>& /*local_M_data*/,
    std::vector<double>& /*local_K_data*/,
    std::vector<double>& /*local_b_data*/,
    std::vector<double>& /*local_Jac_data*/,
    LocalCoupledSolutions const& /*local_coupled_solutions*/)
{
    OGS_FATAL(
        "The assembleWithJacobianForStaggeredScheme() function is not "
        "implemented in"
        " the local assembler.");
}

void LocalAssemblerInterface::computeSecondaryVariable(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, double const t,
    GlobalVector const& x, CoupledSolutionsForStaggeredScheme const* coupled_xs)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);

    if (coupled_xs != nullptr)
        return;

    auto const local_x = x.get(indices);
    computeSecondaryVariableConcrete(t, local_x);
}

void LocalAssemblerInterface::preTimestep(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    double const t, double const delta_t)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    preTimestepConcrete(local_x, t, delta_t);
}

void LocalAssemblerInterface::postTimestep(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    GlobalVector const& x)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    postTimestepConcrete(local_x);
}

void LocalAssemblerInterface::postNonLinearSolver(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    double const t, bool const use_monolithic_scheme)
{
    auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
    auto const local_x = x.get(indices);

    postNonLinearSolverConcrete(local_x, t, use_monolithic_scheme);
}

void LocalAssemblerInterface::writeIntegrationPointDataBinaryInfo(
    std::size_t const mesh_item_id, int const local_vector_size,
    unsigned const integration_point_number, std::ofstream& out)
{
    out.write((char*)(&mesh_item_id), sizeof(std::size_t));
    out.write((char*)(&local_vector_size), sizeof(int));
    out.write((char*)(&integration_point_number), sizeof(unsigned));
}

void LocalAssemblerInterface::checkIntegrationPointDataBinary(
    std::size_t const mesh_item_id, int const local_vector_size,
    unsigned const integration_point_number, std::ifstream& in)
{
    std::size_t original_mesh_item_id = 0;
    in.read((char*)(&original_mesh_item_id), sizeof(std::size_t));
    if (mesh_item_id != original_mesh_item_id)
    {
        OGS_FATAL(
            "The element ID read from the binary file does not match that is "
            "in the restarted computation. Please make sure for a restarted "
            "computation that:\n"
            "\t1. the element ID %d is not changed from the preceding "
            "computation.\n"
            "\t2. or the mesh file is not changed from the preceding "
            "computation.",
            mesh_item_id);
    }

    int vector_size = 0;
    in.read((char*)(&vector_size), sizeof(int));
    if (vector_size != local_vector_size)
    {
        OGS_FATAL(
            "The size of local vector does not match the original one in the "
            "binary data. Reading of the binary file  of the integration point "
            "data stops at element %d",
            mesh_item_id);
    }

    unsigned n_ip = 0;
    in.read((char*)(&n_ip), sizeof(unsigned));
    if (n_ip != integration_point_number)
    {
        OGS_FATAL(
            "The number of integration points does not match the one read from "
            "the binary file of the integration point data for a restarted "
            "computation. Please make sure for a restarted computation that:\n"
            "\t1. the order of numerical integration is not touched from "
            "the preceding computation.\n"
            "\t2. the mesh file is not changed from the preceding "
            "computation.\n Reaing stops at element %d",
            mesh_item_id);
    }
}

}  // namespace ProcessLib
