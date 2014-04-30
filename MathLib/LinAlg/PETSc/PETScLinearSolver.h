/*!
   \file  PETScLinearSolver.h
   \brief Declaration of class PETScLinearSolver, which defines a solver object
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

#ifndef PETSCLINEARSOLVER_H_
#define PETSCLINEARSOLVER_H_

#include "PETScMatrix.h"
#include "PETScVector.h"
#include "PETScLinearSolverOption.h"

namespace MathLib
{

/*!
     A class of linear solver based on PETSc rountines.

*/
class PETScLinearSolver
{
    public:

        /*!
            Constructor.
            \param A           Matrix, cannot be constant.
            \param solver_name Solver name.
            \param pc_name Preconditioner name.
        */
        PETScLinearSolver(PETScMatrix &A, const std::string solver_name, const std::string pc_name);

        ~PETScLinearSolver()
        {
            PCDestroy(_pc);
            KSPDestroy(_solver);
        }

        /*!
            Configure the solver. Call it for re-configuration.
            \param opt Information to configure a solver.
        */

        void setOption(const PETScLinearSolverOption &opt);
        /*!
            Solve a system of equations.
            \param b The right hand of the equations.
            \param x The solutions to be solved.
        */
        void solve(const PETScVector &b, PETScVector &x);

    private:
        KSP *_solver; ///< Slover type.
        PC *_pc;      ///< Preconditioner type.

        /*!
            Set option for preconditioner
            \param set_pc_option Function to set PC options.
            \param opt           Options for solvers and PCs.
        */
        template<typename T_FUNC> void setPC_Option(T_FUNC set_pc_option, const PETScLinearSolverOption &opt);
};

template<typename T_FUNC> void PETScLinearSolver::
setPC_Option(T_FUNC set_pc_option, const PETScLinearSolverOption &opt)
{
    set_pc_option(*_pc, opt);
}

} // end namespace
#endif

