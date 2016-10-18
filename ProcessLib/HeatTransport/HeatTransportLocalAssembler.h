/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   HeatTransportLocalAssembler.h
 *
 */

#ifndef OGS_HEAT_TRANSPORT_LOCALASSEMBLER_H
#define OGS_HEAT_TRANSPORT_LOCALASSEMBLER_H

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HeatTransportMaterialProperties.h"

namespace ProcessLib
{
namespace HeatTransport
{
const unsigned NUM_NODAL_DOF = 1;

class HeatTransportLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtHeatFluxX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtHeatFluxZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class HeatTransportLocalAssembler : public HeatTransportLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

public:
    HeatTransportLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportMaterialProperties& material_propertries)
        : _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _material_properties(material_propertries)
    {
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtHeatFluxX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 0);
        return _heat_fluxes[0];
    }

    std::vector<double> const& getIntPtHeatFluxY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 1);
        return _heat_fluxes[1];
    }

    std::vector<double> const& getIntPtHeatFluxZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_heat_fluxes.size() > 2);
        return _heat_fluxes[2];
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices> _shape_matrices;

    std::vector<std::vector<double>> _heat_fluxes =
        std::vector<std::vector<double>>(
            GlobalDim,
            std::vector<double>(_integration_method.getNumberOfPoints()));

    HeatTransportMaterialProperties& _material_properties;
    double _pore_pressure = 1.e+5;
};

}  // end of namespace
}  // end of namespace

#include "HeatTransportLocalAssembler-impl.h"

#endif /* OGS_HEAT_TRANSPORT_LOCALASSEMBLER_H */
