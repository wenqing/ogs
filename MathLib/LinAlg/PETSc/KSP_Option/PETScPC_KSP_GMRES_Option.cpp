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

#include <logog/include/logog.hpp>

#include "PETScLinearSolverOption.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScPC_KSP_Richards_Option::
PETScPC_KSP_Richards_Option(const boost::property_tree::ptree &option)
    : damping_factor_richards(1.0)
{
    auto damping_factor = option.get_optional<double>("damping_factor");
    damping_factor_richards = *damping_factor;
}

PETScPC_KSP_Chebyshev_Option::
PETScPC_KSP_Chebyshev_Option(const boost::property_tree::ptree &option)
    : emin_chebyshev(0.01), emax_chebyshev(100.0)
{
    auto val = option.get_optional<double>("smallest_eignvalue");
    emin_chebyshev = *val;

    val = option.get_optional<double>("maximum_eignvalue");
    emax_chebyshev = *val;
}

PETScPC_KSP_GMRES_Option::
PETScPC_KSP_GMRES_Option(const boost::property_tree::ptree &option)
    : restart_number_gmres(30), is_modified_gram_schmidt_gmres(false),
      refine_type_gmres(KSP_GMRES_CGS_REFINE_NEVER)
{
    auto val = option.get_optional<double>("restart_number");
    restart_number_gmres = *val;

    boost::optional<bool> bool_vals = option.get_optional<bool>("is_modified_ram_schmidt_orthog");
    is_modified_gram_schmidt_gmres = *bool_vals;

    auto refine_type = option.get_optional<int>("refine_type");
    switch(*refine_type)
    {
        case 0:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_NEVER;
            break;
        case 1:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_IFNEEDED;
            break;
        case 2:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_ALWAYS;
            break;
        default:
            refine_type_gmres = KSP_GMRES_CGS_REFINE_NEVER;
            break;
    }
}

/// Set Chebyshev option
void PETScPC_KSP_GMRES_Option::setOption(KSP &solver)
{
    KSPGMRESSetRestart(solver, restart_number_gmres);

    if(is_modified_gram_schmidt_gmres)
    {
        KSPGMRESSetOrthogonalization(solver, KSPGMRESClassicalGramSchmidtOrthogonalization);
    }

    KSPGMRESSetCGSRefinementType(solver, refine_type_gmres);
}

PETScPC_ILU_Option::
PETScPC_ILU_Option(const boost::property_tree::ptree &option)
    : levels(PETSC_DECIDE), reuse_ordering(false),
      reuse_fill(false), use_in_place(false), allow_diagonal_fill(false)
{
    auto val = option.get_optional<double>("levels");
    levels	= *val;

    auto reuse_order = option.get_optional<bool>("reuse_ordering");
    reuse_ordering = *reuse_order;

    auto rfill = option.get_optional<bool>("reuse_fill");
    reuse_fill = *rfill;

    auto inplane = option.get_optional<bool>("use_in_place");
    use_in_place = *inplane;

    auto low_diag_fill = option.get_optional<bool>("allow_diagonal_fill");
    allow_diagonal_fill = *low_diag_fill;
}

void PETScPC_ILU_Option::setOption(PC &pc)
{
    PCFactorSetLevels(pc, levels);

    if(reuse_ordering)
    {
        PCFactorSetReuseOrdering(pc, PETSC_TRUE);
    }

    if(reuse_fill)
    {
        PCFactorSetReuseFill(pc, PETSC_TRUE);
    }

    if(use_in_place)
    {
        PCFactorSetUseInPlace(pc);
    }

    if(allow_diagonal_fill)
    {
        PCFactorSetAllowDiagonalFill(pc);
    }
}

PETScPC_SOR_Option::
PETScPC_SOR_Option(const boost::property_tree::ptree &option)
    : omega(1.), its(PETSC_DEFAULT),
      lits(PETSC_DEFAULT), type(SOR_FORWARD_SWEEP)
{
    auto val = option.get_optional<double>("omega");
    omega	= *val;

    auto n_its = option.get_optional<int>("local_iterations");
    lits = *n_its;

    n_its = option.get_optional<int>("parallel_iterations");
    its = *n_its;

    auto type_name = option.get_optional<std::string>("type");
    if(type_name->find("pc_sor_symmetric") != std::string::npos)
    {
        type = SOR_FORWARD_SWEEP;
    }
    if(type_name->find("pc_sor_backward") != std::string::npos)
    {
        type = SOR_BACKWARD_SWEEP;
    }
    if(type_name->find("pc_sor_local_forward") != std::string::npos)
    {
        type = SOR_LOCAL_FORWARD_SWEEP;
    }
    if(type_name->find("pc_sor_local_symmetric") != std::string::npos)
    {
        type = SOR_LOCAL_SYMMETRIC_SWEEP;
    }
    if(type_name->find("pc_sor_local_backward") != std::string::npos)
    {
        type = SOR_LOCAL_BACKWARD_SWEEP;
    }
}

void PETScPC_SOR_Option::setOption(PC &pc)
{
    PCSORSetOmega(pc, omega);
    PCSORSetIterations(pc, its, lits);
    PCSORSetSymmetric(pc, type);
}

PETScPC_LU_Option::
PETScPC_LU_Option(const boost::property_tree::ptree &option)
    : mat_type(MATORDERINGNATURAL)
{
    auto type_name = option.get_optional<std::string>("type");
    if(type_name->find("natural") != std::string::npos)
    {
        mat_type = MATORDERINGNATURAL;
    }
    if(type_name->find("nd") != std::string::npos)
    {
        mat_type = MATORDERINGND;
    }
    if(type_name->find("1wd") != std::string::npos)
    {
        mat_type = MATORDERING1WD;
    }
    if(type_name->find("rcm") != std::string::npos)
    {
        mat_type = MATORDERINGRCM;
    }
    if(type_name->find("qmd") != std::string::npos)
    {
        mat_type = MATORDERINGQMD;
    }
    if(type_name->find("rowlength") != std::string::npos)
    {
        mat_type = MATORDERINGROWLENGTH;
    }
}

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

