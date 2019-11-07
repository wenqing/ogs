/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/PhaseFieldExtension.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HydroMechanicalPhaseFieldProcessData.h"
#include "LocalAssemblerInterface.h"

#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace HydroMechanicalPhaseField
{
template <typename BMatricesType, typename ShapeMatrixType, int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;

    typename BMatricesType::KelvinVectorType eps, eps_prev;

    typename BMatricesType::KelvinVectorType sigma_tensile, sigma_compressive,
        sigma_eff;
    double strain_energy_tensile, elastic_energy, pressure, pressure_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C_tensile, C_compressive;
    double integration_weight;

    void pushBackState()
    {
        eps_prev = eps;
        pressure_prev = pressure;
        material_state_variables->pushBackState();
    }

    using Invariants = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value>;

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const /*dt*/,
                                    DisplacementVectorType const& /*u*/,
                                    double const degradation, int split)
    {
        if (split == 0)
        {
            static_cast<
                MaterialLib::Solids::PhaseFieldExtension<DisplacementDim> const&>(
                solid_material)
                .calculateIsotropicDegradedStress(
                    t, x_position, eps, strain_energy_tensile, sigma_tensile,
                    sigma_compressive, C_tensile, C_compressive, sigma_eff,
                    degradation, elastic_energy);
        }
        else if (split == 1)
        {
            static_cast<
                MaterialLib::Solids::PhaseFieldExtension<DisplacementDim> const&>(
                solid_material)
                .calculateDegradedStress(
                    t, x_position, eps, strain_energy_tensile, sigma_tensile,
                    sigma_compressive, C_tensile, C_compressive, sigma_eff,
                    degradation, elastic_energy);
        }
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class HydroMechanicalPhaseFieldLocalAssembler
    : public HydroMechanicalPhaseFieldLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

    HydroMechanicalPhaseFieldLocalAssembler(
        HydroMechanicalPhaseFieldLocalAssembler const&) = delete;
    HydroMechanicalPhaseFieldLocalAssembler(
        HydroMechanicalPhaseFieldLocalAssembler&&) = delete;

    HydroMechanicalPhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicalPhaseFieldProcessData<DisplacementDim>& process_data,
        int const mechanics_related_process_id,
        int const phase_field_process_id,
        int const hydro_process_id)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric),
          _mechanics_related_process_id(mechanics_related_process_id),
          _phase_field_process_id(phase_field_process_id),
          _hydro_process_id(hydro_process_id)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

            static const int kelvin_vector_size =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.C_tensile.setZero(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_compressive.setZero(kelvin_vector_size,
                                          kelvin_vector_size);
            ip_data.sigma_tensile.setZero(kelvin_vector_size);
            ip_data.sigma_compressive.setZero(kelvin_vector_size);
            ip_data.sigma_eff.setZero(kelvin_vector_size);
            ip_data.strain_energy_tensile = 0.0;
            ip_data.elastic_energy = 0.0;
            ip_data.pressure = 0.0;
            ip_data.pressure_prev = 0.0;

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "HydroMechanicalPhaseFieldLocalAssembler: assembly without "
            "Jacobian is not implemented.");
    }

    void assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions) override;

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
        _process_data.width_prev = _process_data.width;
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    void computeFractureNormal(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        CoupledSolutionsForStaggeredScheme const* const cpl_xs) override;

    void computeFractureWidth(
        std::size_t mesh_item_id,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        double const t, CoupledSolutionsForStaggeredScheme const* const cpl_xs,
        MeshLib::Mesh const& mesh) override;

    void findNeighborElement(MeshLib::Element const& current_ele,
                             GeoLib::LineSegment& LIntegral,
                             MeshLib::Element const*& neighbor_ele,
                             Eigen::Vector3d& intersectionMidPoint,
                             std::size_t last_visited) override;

/*    void findHostElement(MeshLib::Element const& current_ele,
                         Eigen::Vector3d pnt_end,
                         MeshLib::Element const*& neighbor_ele) override;*/

    /*    void computeEnergy(
            std::size_t mesh_item_id,
            std::vector<
                std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
                dof_tables,
            GlobalVector const& x, double const t, double& elastic_energy,
            double& surface_energy, double& pressure_work,
            bool const use_monolithic_scheme,
            CoupledSolutionsForStaggeredScheme const* const cpl_xs) override;*/

private:
    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma_eff = _ip_data[ip].sigma_eff;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_eff);
        }

        return cache;
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
        }

        return cache;
    }

    void assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    void assembleWithJacobianForHydroProcessEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    void assembleWithJacobianForPhaseFieldEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    HydroMechanicalPhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunction::NPOINTS;
    static const int phasefield_index = ShapeFunction::NPOINTS;
    static const int phasefield_size = ShapeFunction::NPOINTS;
    static const int displacement_index = 2 * ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;

    /// ID of the processes that contains mechanical process.
    int const _mechanics_related_process_id;

    /// ID of phase field process.
    int const _phase_field_process_id;

    /// ID of hydro process.
    int const _hydro_process_id;
};

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib

#include "HydroMechanicalPhaseFieldFEM-impl.h"
