/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>

#include "MeshLib/Elements/Hex.h"

namespace NumLib
{
/**
 *  Shape function for a 20-nodes hex element in natural coordinates
 *
 */
class ShapeHex20
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  rst natural coordinates (r,s,t)
     * @param [out] N   a vector of calculated shape functions
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X& rst, T_N& N);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  rst natural coordinates (r,s,t)
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X& rst, T_N& dN);

    static constexpr std::array reference_element_centre = {0.0, 0.0, 0.0};

    using MeshElement = MeshLib::Hex20;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
    static constexpr int ORDER = 2;
};

}  // namespace NumLib

#include "ShapeHex20-impl.h"
