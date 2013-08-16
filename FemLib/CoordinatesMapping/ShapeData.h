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


#ifndef SHAPEDATA_H_
#define SHAPEDATA_H_

#include "MathLib/DataType.h"

namespace FemLib
{

/**
 * \brief Coordinate mapping matrices at particular location
 *
 * Mapping is done between physical coordinates (x,y,z) and natural coordinates (r,s,t)
 */
struct ShapeData
{
    /// shape function N(r)
    MathLib::LocalRowVector N;

    /// gradient of shape functions, dN(r)/dr
    MathLib::LocalMatrix dNdr;

    /// gradient of shape functions, dN(r)/dx
    MathLib::LocalMatrix dNdx;

    /// Jacobian matrix, J=dx/dr
    MathLib::LocalMatrix J;

    /// inverse of the Jacobian
    MathLib::LocalMatrix invJ;

    /// determinant of the Jacobian
    double detJ;

    /**
     *
     * @param dim       Dimension of the physical coordinates
     * @param n_nodes   The number of element nodes
     */
    ShapeData(std::size_t dim, std::size_t n_nodes)
    : N(MathLib::LocalRowVector::Zero(n_nodes)), dNdr(MathLib::LocalMatrix::Zero(dim, n_nodes)),
      dNdx(MathLib::LocalMatrix::Zero(dim, n_nodes)), J(MathLib::LocalMatrix::Zero(dim, dim)),
      invJ(MathLib::LocalMatrix::Zero(dim, dim)), detJ(.0)
    {}

    ~ShapeData() {}

    void setZero()
    {
        setZero(N);
        setZero(dNdr);
        setZero(dNdx);
        setZero(J);
        setZero(invJ);
        detJ = .0;
    }

private:
    void setZero(MathLib::LocalMatrix &mat)
    {
        mat.setZero(mat.rows(), mat.cols());
    }

    void setZero(MathLib::LocalRowVector &vec)
    {
        vec.setZero(vec.size());
    }

};

}

#endif //SHAPEDATA_H_
