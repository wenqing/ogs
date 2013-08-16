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


#ifndef NATURALCOORDINATESMAPPING_H_
#define NATURALCOORDINATESMAPPING_H_


#include "MeshLib/Elements/Element.h"

#include "ShapeData.h"


namespace FemLib
{

/**
 * \brief Mapping element shapes to natural coordinates
 *
 * This class is responsible for mapping element shapes in physical coordinates
 * (x,y,z) to those in natural coordinates (r,s,t).
 * - Given physical coordinates should correspond to dimensions of the element, 
 *   i.e. (x,y) for triangles
 * - _nodes_coords(r,s,t) = N(r,s,t) x_i
 */
template <class T_SHAPE_FUNC>
class NaturalCoordinatesMapping
{
public:
    NaturalCoordinatesMapping() : _ele(nullptr) {}

    explicit NaturalCoordinatesMapping(const MeshLib::Element &ele);

    ~NaturalCoordinatesMapping() {}

    /// reset a mesh element
    void reset(const MeshLib::Element &ele);

    /// compute mapping matrices at the given location in natural coordinates
    ///
    /// @param natural_pt
    /// @param shape
    void computeMappingMatrices(const double* natural_pt, ShapeData &shape) const;

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{_nodes_coords} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param shape
    /// @param physical_pt
    void mapToPhysicalCoordinates(const ShapeData &shape, double* physical_pt) const;

    /// compute natural coordinates at the given physical coordinates.
    /// Assuming \f$ r=0 \f$ at \f$ _nodes_coords = \bar{_nodes_coords}_{avg} \f$, natural coordinates can be calculated as
    /// \f[
    ///    \mathbf{r} = (\mathbf{J}^-1)^T * (\mathbf{_nodes_coords} - \bar{\mathbf{_nodes_coords}}_{avg})
    /// \f]
    ///
    /// @param shape
    /// @param physical_pt
    /// @param natural_pt
    void mapToNaturalCoordinates(const ShapeData &shape, const double* physical_pt, double* natural_pt) const;

private:
    const MeshLib::Element* _ele;
    MathLib::LocalMatrix _nodes_coords;
};

} //end namespace

#include "NaturalCoordinatesMapping.tpp"


#endif //NATURALCOORDINATESMAPPING_H_
