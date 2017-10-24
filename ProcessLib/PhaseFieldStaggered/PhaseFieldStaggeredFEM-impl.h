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
#include "NumLib/Function/Interpolation.h"

namespace ProcessLib
{
namespace PhaseFieldStaggered
{

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>

void LocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
     assemblePhaseField(SpatialPosition& pos,
        std::vector<double> const& local_x, std::vector<double> const& local_u,
        std::vector<double> const& strain_energy_tensile_ips, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_rhs_data)
{
     auto const local_matrix_size = local_x.size();
     // This assertion is valid only if all nodal d.o.f. use the same shape
     // matrices.
     assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

     auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
         local_M_data, local_matrix_size, local_matrix_size);

     auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
         local_K_data, local_matrix_size, local_matrix_size);

     auto local_rhs = MathLib::createZeroedVector<RhsVector>(
         local_rhs_data, local_matrix_size);

     // Not sure about getting d
     auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
         phasefield_size> const>(local_x.data(), phasefield_size);

     const auto local_u_vec =
         MathLib::toVector<NodalVectorType>(local_u, local_matrix_size);

     unsigned const n_integration_points =
         _integration_method.getNumberOfPoints();

     assert(strain_energy_tensile_ips.size() == n_integration_points);


     for (unsigned ip = 0; ip < n_integration_points; ip++)
     {
         auto const& w = _integration_method.getWeightedPoint(ip);
         auto const& N = _shape_matrices[ip].N;
         auto const& dNdx = _shape_matrices[ip].dNdx;
         double const d_ip = N.dot(d);

         auto const x_coord =
             interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                 _element, N);

         double const gc = _process_data.crack_resistance(t, x_position)[0];
         double const ls = _process_data.crack_length_scale(t, x_position)[0];

         // phasefield equation.
         local_K.noalias() += (gc * (0.5 * N.transpose() * N / ls + 2 * dNdx.transpose() * dNdx * ls )
                           + N.transpose() * N * strain_energy_tensile_ips ) * w;

         local_rhs.noalias() -= (2 * dNdx.transpose() * dNdx * ls* d +
                            N.transpose() * d_ip * 2 * strain_energy_tensile_ips -
                            N.transpose() * 0.5 * gc / ls * (1 - d_ip)) * w;
     }

}

void PhaseFieldStaggeredLocalAssemblerData <typename ShapeFunction, typename IntegrationMethod,
    unsigned GlobalDim>::assembleWithCoupledTerm(
    double const t, std::vector<double> const& local_x,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& /*local_b_data*/,
    LocalCouplingTerm const& coupled_term)
{

    for (auto const& coupled_process_pair : coupled_term.coupled_processes)
    {
        if (coupled_process_pair.first ==
                std::type_index(typeid(ProcessLib::PhaseFieldSmallDeformation)))
        {
            assert(dynamic_cast<ProcessLib::PhaseFieldSmallDeformation::PhaseFieldSmallDeformationProcess*>
                   (&(coupled_process_pair.second)) != nullptr);

            auto const& pcs =
                    static_cast<const ProcessLib::PhaseFieldSmallDeformation::PhaseFieldSmallDeformationProcess const&>
                    (coupled_process_pair.second);

            auto const strain_energy_tensile_ips = pcs.getIntStrainEnergyTensile(_element.getID());

            const auto local_u =
                    coupled_term.local_coupled_xs.at(coupled_process_pair.first);
            SpatialPosition pos;
            pos.setElementID(_element.getID());

            assemblePhaseField(pos, local_x, local_u, strain_energy_tensile_ips, local_M_data, local_K_data,
               local_rhs_data)
        }

    }

}


}  // namespace PhaseFieldStaggered
}  // namespace ProcessLib
