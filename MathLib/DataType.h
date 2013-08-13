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

#ifndef MATHLIB_DATATYPE_H_
#define MATHLIB_DATATYPE_H_

#include <Eigen>

namespace MathLib
{

/// Local dense matrix type (row-majored)
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;

/// Local dense vector type
typedef Eigen::VectorXd LocalVector;

typedef Eigen::RowVectorXd LocalRowVector;

}

#endif // MATHLIB_DATATYPE_H_
