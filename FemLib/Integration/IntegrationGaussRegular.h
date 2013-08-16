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

#ifndef INTEGRATIONGAUSSREGULAR_H_
#define INTEGRATIONGAUSSREGULAR_H_

#include <tuple>

#include "MathLib/Integration/GaussLegendre.h"

#include "IIntegrationGauss.h"

namespace FemLib
{

/**
 * \brief Gauss quadrature rule for quad
 *
 */
template <std::size_t N_DIM>
class IntegrationGaussRegular : public IIntegrationGauss
{
public:
    /**
     * Construct this object with the given sampling level
     *
     * @param n_sampl_level     the sampling level (default 2)
     */
    explicit IntegrationGaussRegular(std::size_t n_sampl_level = 2)
    : _n_sampl_level(n_sampl_level), _n_sampl_pt(0)
    {
        this->setSamplingLevel(n_sampl_level);
    }

    /**
     * Change the sampling level
     *
     * @param n_sampl_level     the sampling level
     */
    virtual void setSamplingLevel(std::size_t n_sampl_level)
    {
        this->_n_sampl_pt = std::pow(n_sampl_level, N_DIM);
        this->_n_sampl_level = n_sampl_level;
    }

    /// return current sampling level
    virtual std::size_t getSamplingLevel() const {return _n_sampl_level;}

    /// return the number of sampling points
    virtual std::size_t getNPoints() const {return _n_sampl_pt;}

    /**
     * get a sampling point
     *
     * @param igp   the sampling point index
     * @param x     coordinates of the sampling point
     */
    virtual void getPoint(std::size_t igp, double* x) const
    {
        getPoint(getSamplingLevel(), igp, x);
    }

    /// return weight at the given sampling point
    virtual double getWeight(std::size_t igp) const
    {
        return getWeight(getSamplingLevel(), igp);
    }

    /// get position indexes
    static std::tuple<std::size_t, std::size_t, std::size_t> getPosition(std::size_t /*nGauss*/, std::size_t /*igp*/)
    {
        return std::make_tuple(0u, 0u, 0u);
    }

    /// get a sampling point
    static void getPoint(std::size_t nGauss, std::size_t igp, double* x);

    /// return weight at the given sampling point
    static double getWeight(std::size_t nGauss, std::size_t igp);

private:
    std::size_t _n_sampl_level;
    std::size_t _n_sampl_pt;
};

} //namespace

#include "IntegrationGaussRegular.tpp"

#endif //INTEGRATIONGAUSSREGULAR_H_
