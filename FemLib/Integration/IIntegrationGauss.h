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


#ifndef IINTEGRATIONGAUSS_H_
#define IINTEGRATIONGAUSS_H_

#include <cstddef>

namespace FemLib
{

/**
 * \brief interface of Gauss quadrature methods
 */
class IIntegrationGauss
{
public:
    virtual ~IIntegrationGauss() {}

    /**
     * Set the sampling level
     *
     * @param n_sampl_level
     */
    virtual void setSamplingLevel(std::size_t /*n_sampl_level*/) {};

    /// return current sampling level
    virtual std::size_t getSamplingLevel() const {return 0;}

    /// return the total number of sampling points
    virtual std::size_t getNPoints() const = 0;

    /**
     * get a sampling point coordinates
     *
     * @param pt_id    the sampling point id
     * @param x        coordinates
     */
    virtual void getPoint(std::size_t pt_id, double* x) const = 0;

    /**
     * get a weight of the sampling point
     *
     * @param pt_id    the sampling point id
     * @return weight
     */
    virtual double getWeight(std::size_t pt_id) const = 0;
};

} // end namespace

#endif //IINTEGRATIONGAUSS_H_
