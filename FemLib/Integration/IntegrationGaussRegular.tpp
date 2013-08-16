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


#include "MathLib/Integration/GaussLegendre.h"

namespace FemLib
{

template <>
inline void IntegrationGaussRegular<1>::getPoint(std::size_t nGauss, std::size_t igp, double* x)
{
    x[0] = MathLib::GaussLegendre::getPoint(nGauss, igp);
}

template <>
inline double IntegrationGaussRegular<1>::getWeight(std::size_t nGauss, std::size_t igp)
{
    return  MathLib::GaussLegendre::getWeight(nGauss, igp);
}

template <>
inline std::tuple<std::size_t, std::size_t, std::size_t>
IntegrationGaussRegular<2>::getPosition(std::size_t nGauss, std::size_t igp)
{
    std::size_t gp_r = igp / nGauss;
    std::size_t gp_s = igp % nGauss;
    return std::make_tuple(gp_r, gp_s, 0u);
}

template <>
inline void IntegrationGaussRegular<2>::getPoint(std::size_t nGauss, std::size_t igp, double* x)
{
    auto pos = getPosition(nGauss, igp);
    x[0] = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    x[1] = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
}

template <>
inline double IntegrationGaussRegular<2>::getWeight(std::size_t nGauss, std::size_t igp)
{
    auto pos = getPosition(nGauss, igp);
    return  MathLib::GaussLegendre::getWeight(nGauss, std::get<0>(pos))
            * MathLib::GaussLegendre::getWeight(nGauss, std::get<1>(pos));
}

template <>
inline std::tuple<std::size_t, std::size_t, std::size_t>
IntegrationGaussRegular<3>::getPosition(std::size_t nGauss, std::size_t igp)
{
    std::size_t gp_r = igp / (nGauss * nGauss);
    std::size_t gp_s = igp % (nGauss * nGauss);
    std::size_t gp_t = gp_s % nGauss;
    gp_s /= nGauss;
    return std::make_tuple(gp_r, gp_s, gp_t);
}

template <>
inline void IntegrationGaussRegular<3>::getPoint(std::size_t nGauss, std::size_t igp, double* x)
{
    auto pos = getPosition(nGauss, igp);
    x[0] = MathLib::GaussLegendre::getPoint(nGauss, std::get<0>(pos));
    x[1] = MathLib::GaussLegendre::getPoint(nGauss, std::get<1>(pos));
    x[2] = MathLib::GaussLegendre::getPoint(nGauss, std::get<2>(pos));
}

template <>
inline double IntegrationGaussRegular<3>::getWeight(std::size_t nGauss, std::size_t igp)
{
    auto pos = getPosition(nGauss, igp);
    return  MathLib::GaussLegendre::getWeight(nGauss, std::get<0>(pos))
            * MathLib::GaussLegendre::getWeight(nGauss, std::get<1>(pos))
            * MathLib::GaussLegendre::getWeight(nGauss, std::get<2>(pos));
}

} //namespace

