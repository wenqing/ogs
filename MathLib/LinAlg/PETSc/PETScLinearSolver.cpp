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
using boost::property_tree::ptree;

PETScLinearSolver::PETScLinearSolver(PETScMatrix &A,
                                     const boost::property_tree::ptree &option)
{
    KSPCreate(PETSC_COMM_WORLD, &_solver);
    KSPSetOperators(_solver, A.getRawMatrix(), A.getRawMatrix(), DIFFERENT_NONZERO_PATTERN);

    boost::optional<ptree> pt_solver = option.get_child("linear_solver");
    if(!pt_solver)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n*** PETSc linear solver is not specified, bcgs + bjacobi are used.\n");
        return;
    }

    // Preconditioners:
    boost::optional<ptree> pt_pc = option.get_child("preconditioner");
    if(!pt_pc)
    {
        PetscPrintf(PETSC_COMM_WORLD,"\n*** PETSc preconditioner is not specified, bjacobi is used.");
        return;
    }

    // Base configuration
    PETScLinearSolverOption opt(*pt_solver, *pt_pc);
          
    opt.setOption(_solver, _pc);

// To be extend:

    //----------------------------------------------------------------------
    // Specific configuration, solver
    boost::optional<ptree> pt_solver_spec = pt_solver->get_child("Richards");
    if(pt_solver_spec)
    {
        PETScPC_KSP_Richards_Option ksp_opt(*pt_solver_spec);
        setKSP_Option(ksp_opt);
    }

    pt_solver_spec = pt_solver->get_child("Chebyshev");
    if(pt_solver_spec)
    {
        PETScPC_KSP_Chebyshev_Option ksp_opt(*pt_solver_spec);
        setKSP_Option(ksp_opt);
    }

    pt_solver_spec = pt_solver->get_child("gmres");
    if(pt_solver_spec)
    {
        PETScPC_KSP_GMRES_Option ksp_opt(*pt_solver_spec);
        setKSP_Option(ksp_opt);
    }

    //----------------------------------------------------------------------
    // Specific configuration, preconditioner
    // ILU or ICC
    boost::optional<ptree> pt_pc_spec = pt_pc->get_child("ilu");
    if(pt_pc_spec)
    {
        PETScPC_ILU_Option pc_opt(*pt_pc_spec);
        setPC_Option(pc_opt);
    }
    else
    {
        pt_pc_spec = pt_pc->get_child("icc");
        if(pt_pc_spec)
        {
            PETScPC_ILU_Option pc_opt(*pt_pc_spec);
            setPC_Option(pc_opt);
        }
    }

    pt_pc_spec = pt_pc->get_child("sor");
    if(pt_pc_spec)
    {
        PETScPC_SOR_Option pc_opt(*pt_pc_spec);
        setPC_Option(pc_opt);
    }

    pt_pc_spec = pt_pc->get_child("lu");
    if(pt_pc_spec)
    {
        PETScPC_LU_Option pc_opt(*pt_pc_spec);
        setPC_Option(pc_opt);
    }

    pt_pc_spec = pt_pc->get_child("asm");
    if(pt_pc_spec)
    {
        PETScPC_ASM_Option pc_opt(*pt_pc_spec);
        setPC_Option(pc_opt);
    }

    pt_pc_spec = pt_pc->get_child("amg");
    if(pt_pc_spec)
    {
        PETScPC_AMG_Option pc_opt(*pt_pc_spec);
        setPC_Option(pc_opt);
    }  
           
    //
    KSPSetFromOptions(_solver);  // set running time option
}

void PETScLinearSolver::solve(const PETScVector &b, PETScVector &x)
{
// define TEST_MEM_PETSC
#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    KSPSolve(_solver, b.getData(), x.getData());

    KSPConvergedReason reason;
    KSPGetConvergedReason(_solver, &reason);

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
        KSPGetType(_solver, &slv_type);
        PCGetType(_pc, &prc_type);

        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
        PetscPrintf(PETSC_COMM_WORLD, "\nLinear solver %s with %s preconditioner",
                    slv_type, prc_type);

        int its;
        KSPGetIterationNumber(_solver, &its);
        PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n", its);
        PetscPrintf(PETSC_COMM_WORLD,"\n================================================");
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");

#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

} //end of namespace

