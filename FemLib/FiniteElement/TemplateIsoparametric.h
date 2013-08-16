/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef TEMPLATEISOPARAMETRIC_H_
#define TEMPLATEISOPARAMETRIC_H_


#include <vector>
#include <cassert>

#include "FemLib/CoordinatesMapping/NaturalCoordinatesMapping.h"

namespace FemLib
{

/**
 * \brief Template class for isoparametric elements
 *
 * \tparam T_ELEMENT        Mesh element class
 * \tparam T_SHAPE          Shape function
 * \tparam T_INTEGRAL       Integration method
 * \tparam T_EXTRAPOLATE    Extrapolation method
 */
template <
    class T_MESH_ELEMENT,
    class T_SHAPE,
    class T_INTEGRAL,
    class T_EXTRAPOLATE
    >
class TemplateIsoparametric
{
public:
    /**
     * Constructor without specifying a mesh element. resetMeshElement() must be called afterwards.
     *
     * @param integration_points_level      the sampling level (default 2)
     */
    explicit TemplateIsoparametric(std::size_t integration_points_level=2)
    : _ele(nullptr)
    {
        this->_integration.setSamplingLevel(integration_points_level);
    }

    /**
     * Construct this object for the given mesh element.
     *
     * @param e                             Mesh element object
     * @param integration_points_level      the sampling level (default 2)
     */
    TemplateIsoparametric(const T_MESH_ELEMENT &e, std::size_t integration_points_level=2)
    : _ele(&e)
    {
        this->resetMeshElement(e);
        this->_integration.setSamplingLevel(integration_points_level);
    }

    ///
    ~TemplateIsoparametric() {}

    /// return current mesh element
    const T_MESH_ELEMENT* getMeshElement() const {return _ele;}

    /// reset a mesh element
    void resetMeshElement(const T_MESH_ELEMENT &e)
    {
        this->_ele = &e;
        this->_mapping.reset(e);
    }

    /// return an integration method
    const T_INTEGRAL& getIntegrationMethod() const {return _integration;}

    /**
     * compute shape functions
     *
     * @param x         point in natural coordinates
     * @param shape     evaluated shape data
     */
    void computeShapeFunctions(const double *x, ShapeData &shape) const
    {
        _mapping.computeMappingMatrices(x, shape);
    }

    /// make interpolation from nodal values
    double interpolate(const ShapeData& shape, MathLib::LocalVector &nodal_values) const
    {
        assert(nodal_values.size()==this->_ele->getNNodes());
        return shape.N.dot(nodal_values);
    }

    /// extrapolate integration point values to nodal values
    void extrapolate(const MathLib::LocalVector &gp_values, MathLib::LocalVector &nodal_values) const
    {
        T_EXTRAPOLATE::extrapolate(gp_values, nodal_values);
    }

private:
    const T_MESH_ELEMENT* _ele;
    NaturalCoordinatesMapping<T_SHAPE> _mapping;
    T_INTEGRAL _integration;
};

} //end

#endif //TEMPLATEISOPARAMETRIC_H_
