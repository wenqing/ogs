/*!
   \file  PETScLinearSolver.cpp
   \brief Definition of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScLinearSolver.h"

namespace MathLib
{
// Set ILU/ICC preconditioner options
void setPETSC_PC_OptionILU(PC pc, const PETScLinearSolverOption &opt);

// Set SOR/SSOR preconditioner options
void setPETSC_PC_OptionSOR(PC pc, const PETScLinearSolverOption &opt);

PETScLinearSolver::PETScLinearSolver(PETScMatrix &A,
                                     const std::string solver_name, const std::string pc_name)
    : _solver(NULL), _pc(NULL)
{
    KSPCreate(PETSC_COMM_WORLD, _solver);
    KSPSetOperators(*_solver, A.getRawMatrix(), A.getRawMatrix(), DIFFERENT_NONZERO_PATTERN);

    KSPSetType(*_solver, solver_name.c_str());
    KSPGetPC(*_solver, _pc);
    PCSetType(*_pc, pc_name.c_str());
}

void PETScLinearSolver::setOption(const PETScLinearSolverOption &opt)
{
    // --------------------------------------------------------------
    // Preconditioner:
    if(_pc)
    {
        PCDestroy(_pc);
    }
    KSPSetType(*_solver, opt.solver_name.c_str());
    KSPGetPC(*_solver, _pc);
    PCSetType(*_pc, opt.pc_name.c_str());
    KSPSetPCSide(*_solver, opt.preco_side);

    if(   opt.pc_name.find("ilu") != std::string::npos
            || opt.pc_name.find("icc") != std::string::npos )
    {
        setPC_Option(setPETSC_PC_OptionILU, opt);
    }

    if(opt.pc_name.find("sor") != std::string::npos)
    {
        setPC_Option(setPETSC_PC_OptionSOR, opt);
    }


    // --------------------------------------------------------------
    // Solver:
    if(opt.solver_name.find("richardson") != std::string::npos)
    {
        KSPRichardsonSetScale(*_solver, opt.damping_factor_richards);
    }

    if(opt.solver_name.find("chebychev") != std::string::npos)
    {
        KSPChebyshevSetEigenvalues(*_solver, opt.emax_chebyshev, opt.emin_chebyshev);
    }

    if(opt.solver_name.find("gmres") != std::string::npos)
    {
        KSPGMRESSetRestart(*_solver, opt.restart_number_gmres);

        if(opt.is_modified_gram_schmidt_gmres)
        {
            KSPGMRESSetOrthogonalization(*_solver, KSPGMRESClassicalGramSchmidtOrthogonalization);
        }

        KSPGMRESSetCGSRefinementType(*_solver, opt.refine_type_gmres);
    }

    KSPSetTolerances(*_solver, opt.rtol, opt.atol, opt.dtol, opt.max_it);
    KSPSetFromOptions(*_solver);  // set running time option
}


void PETScLinearSolver::solve(const PETScVector &b, PETScVector &x)
{
// define TEST_MEM_PETSC
#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    KSPSolve(*_solver, b.getData(), x.getData());

    KSPConvergedReason reason;
    KSPGetConvergedReason(*_solver, &reason);

    if(reason == KSP_DIVERGED_INDEFINITE_PC)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite preconditioner;\n");
        PetscPrintf(PETSC_COMM_WORLD,"Run the executable again but with -pc_factor_shift_positive_definite option.\n");
    }
    else if(reason < 0)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should not happen.\n");
    }
    else
    {
        const char *slv_type;
        const char *prc_type;
        KSPGetType(*_solver, &slv_type);
        PCGetType(*_pc, &prc_type);

        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
        PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                    slv_type, prc_type);

        int its;
        KSPGetIterationNumber(*_solver, &its);
        PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n", its);
        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");

#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

void setPETSC_PC_OptionILU(PC pc, const PETScLinearSolverOption &opt)
{
    PCFactorSetLevels(pc, opt.pc_ilu.levels);

    if(opt.pc_ilu.reuse_ordering)
    {
        PCFactorSetReuseOrdering(pc, PETSC_TRUE);
    }

    if(opt.pc_ilu.reuse_fill)
    {
        PCFactorSetReuseFill(pc, PETSC_TRUE);
    }

    if(opt.pc_ilu.use_in_place)
    {
        PCFactorSetUseInPlace(pc);
    }

    if(opt.pc_ilu.allow_diagonal_fill)
    {
        PCFactorSetAllowDiagonalFill(pc);
    }
}

void setPETSC_PC_OptionSOR(PC pc, const PETScLinearSolverOption &opt)
{
    PCSORSetOmega(pc, opt.pc_sor.omega);
    PCSORSetIterations(pc, opt.pc_sor.its, opt.pc_sor.lits);
    PCSORSetSymmetric(pc, opt.pc_sor.type);
}

} //end of namespace

