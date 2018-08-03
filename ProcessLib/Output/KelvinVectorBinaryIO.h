/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   KelvinVectorBinaryIO.h
 *
 * Created on August 3, 2018, 2:56 PM
 */
#pragma once

#include <fstream>
#include <Eigen/Dense>

#include "MathLib/KelvinVector.h"

namespace ProcessLib
{
template <typename KelvinVector>
void writeKelvinVectorBinary(std::ofstream& out, const KelvinVector& vector)
{
    const typename KelvinVector::Index size = vector.rows() * vector.cols();
    out.write((char*)(&size), sizeof(typename KelvinVector::Index));
    out.write((char*)vector.data(),
              size * sizeof(typename KelvinVector::Scalar));
}

template <typename KelvinVector>
void readKelvinVectorBinary(std::ifstream& in, KelvinVector& vector)
{
    typename KelvinVector::Index size = 0;
    in.read((char*)(&size), sizeof(typename KelvinVector::Index));
    vector.setZero(size);
    in.read((char*)vector.data(), size * sizeof(typename KelvinVector::Scalar));
}

}  // end of namespace
