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


#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include <vector>
#include <cstddef>

namespace MathLib
{

/**
 * \brief Gauss-Legendre quadrature method
 *
 */
class GaussLegendre
{
public:
    /**
     * return a sampling point
     *
     * \param n_sample_points   the number of sampling points
     * \param point_id          point index
     * \return x
     */
    static double getPoint(std::size_t n_sample_points, std::size_t point_id);
    

    /**
     * return a weigth of the sampling point
     *
     * \param n_sample_points   the number of sampling points
     * \param point_id          point index
     * \return weight
     */
    static double getWeight(std::size_t n_sample_points, std::size_t point_id);

};

}

#endif // GAUSSLEGENDRE_H_

