/*!
   \file  PETScLinearSolverOption.cpp
   \brief Define members of PETScLinearSolverOption

   \author Wenqing Wang
   \date 02-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScLinearSolverOption.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScLinearSolverOption::
PETScLinearSolverOption(const boost::property_tree::ptree &ksp_option,
                        const boost::property_tree::ptree &pc_option)
    : _solver_name("bcgs"), _pc_name("bjacobi"), _preco_side(PC_LEFT),
      _max_it(2000), _rtol(1.e-5), _atol(PETSC_DEFAULT), _dtol(PETSC_DEFAULT)
{
    auto solver_type = ksp_option.get_optional<std::string>("solver_type");
    if(solver_type)
    {
        _solver_name = *solver_type;
    }

    auto max_iteration_step = ksp_option.get_optional<int>("max_iteration_step");
    if(max_iteration_step)
    {
        _max_it = *max_iteration_step;
    }

    auto error_tolerance = ksp_option.get_optional<double>("relative_tolerance");
    if(error_tolerance)
    {
        _rtol = *error_tolerance;
    }

    error_tolerance = ksp_option.get_optional<double>("absolute_tolerance");
    if(error_tolerance)
    {
        _atol = *error_tolerance;
    }

    error_tolerance = ksp_option.get_optional<double>("relative_increase_residual");
    if(error_tolerance)
    {
        _dtol = *error_tolerance;
    }

    // Preconditioners:
    auto pc_side = pc_option.get_optional<std::string>("precond_side");
    if(pc_side)
    {
        if(pc_side->find("right") != std::string::npos)
            _preco_side = PC_RIGHT;
        if(pc_side->find("symmetric") != std::string::npos)
            _preco_side = PC_SYMMETRIC;
    }
}

} // end namespace

