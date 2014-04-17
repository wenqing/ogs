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
            \param A   Matrix.
            \param opt Configuration options for this solver.
        */
        /// Constructor
        PETScLinearSolver(const PETScMatrix &A,
                          const PETScLinearSolverOption opt = PETScLinearSolverOption() );

        ~PETScLinearSolver()
        {
            KSPDestroy(&lsolver);
        }

        /*!
            Configure the solver. Call it for re-configuration.
            \param opt Information to configure a solver.
        */

        void setOption(const PETScLinearSolverOption opt);
        /*!
            Solve a system of equations.
            \param b The right hand of the equations.
            \param x The solutions to be solved.
        */
        void solve(const PETScVector &b, PETScVector &x);

    private:
        KSP _solver; ///< Slover type.
        PC _pc;      ///< Preconditioner type.
};

} // end namespace
#endif

