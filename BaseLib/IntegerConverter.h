/*!
  \file IntegerConverter.h
  \brief  Convert integer types for c style functions like those used in MPI

  \author Wenqing Wang
  \date   2015.01

  \copyright
  Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef INTEGER_CONVERTER_H
#define INTEGER_CONVERTER_H

#include<climits>

namespace BaseLib
{

/// Convert any integer to int type integer
template<typename T_INT> int toInt(const T_INT int_var)
{
    if( int_var > std::numeric_limits<int>::max() )
    {
        throw std::runtime_error("Error: Integer exceeds its maximum limit");
    }

    return static_cast<int>(int_var);
}

}

#endif
