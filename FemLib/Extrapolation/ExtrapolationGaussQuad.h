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


#ifndef EXTRAPOLATIONGAUSSQUAD_H_
#define EXTRAPOLATIONGAUSSQUAD_H_

#include "MathLib/DataType.h"

namespace FemLib
{

/**
 * \brief Extrapolation of Gauss point values to nodal values
 */
class ExtrapolationGaussQuad
{
public:
    /**
     * extrapolate gauss point values to nodal values
     *
     * @param gp_values     a vector of gauss point values
     * @param nodal_values  a vector of extrapolated nodal values
     */
    static void extrapolate(const MathLib::LocalVector &gp_values, MathLib::LocalVector &nodal_values);

private:
    static std::size_t getNodeIndexOfGaussQuad(std::size_t nGaussLevel, std::size_t igp);
    static double calculateXi_p(std::size_t nGaussLevel);
    static void getExtrapolatedPoints(std::size_t nodeIdOfGaussQuad, double Xi_p, double* pt_in_natural);

};

} // end namespace

#endif //EXTRAPOLATIONGAUSSQUAD_H_
