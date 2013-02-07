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


#include "GaussLegendre.h"

namespace MathLib
{

double GaussLegendre::getPoint(std::size_t n_sample_points, std::size_t point_id)
{
    switch (n_sample_points)
    {
    case 1:
        return 0.0;
    case 2:
        switch (point_id)
        {
        case 0:
            return 0.577350269189626;
        case 1:
            return -0.577350269189626;
        }
        break;
    case 3:
        switch (point_id)
        {
        case 0:
            return 0.774596669241483;
        case 1:
            return 0.0;
        case 2:
            return -0.774596669241483;
        }
        break;
    case 4:
        switch (point_id)
        {
        case 0:
            return 0.861136311594053;
        case 1:
            return 0.339981043584856;
        case 2:
            return -0.339981043584856;
        case 3:
            return -0.861136311594053;
        }
        break;
    }
    return 0.0;
}

double GaussLegendre::getWeight(std::size_t n_sample_points, std::size_t point_id)
{
    switch (n_sample_points)
    {
    case 1:
        return 2.0;
    case 2:
        return 1.0;
        break;
    case 3:
        switch (point_id)
        {
        case 0:
        case 2:
            return 0.555555555555556;
        case 1:
            return 0.888888888888889;
        }
        break;
    case 4:
        switch (point_id)
        {
        case 0:
            return 0.347854845137454;
        case 3:
        case 1:
        case 2:
            return 0.652145154862546;
        }
        break;
    }
    return 0.0;
}


} //namespace
