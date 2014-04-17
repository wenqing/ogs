/*!
   \file  PETScLinearSolverOption.h
   \brief Define the configuration data for the PETSc linear solver.

   \author Wenqing Wang
   \date 02-2014

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PETSCLINEARSOLVEROPTION_H_
#define PETSCLINEARSOLVEROPTION_H_

#include <string>

#include <petscksp.h>

namespace MathLib
{
/*!
   \brief This a struct data containing configuration data
          used to configure a procondtioner and a linear solver of PETSc.
*/
struct PETScLinearSolverOption
{
    PETScLinearSolverOption() :  max_it(2000), rtol(1.e-5),
        atol(PETSC_DEFAULT), dtol(PETSC_DEFAULT),
        solver_name("bcgs"), solver_name("bjacobi"),
        damping_factor_richards(1.0),
        emin_chebyshev(0.01), emax_chebyshev(100.0),
        restart_number_gmres(30), is_gram_schmidt_orthog_gmres(true),
        refine_type_gmres(KSP_GMRES_CGS_REFINE_NEVER),
        righ_side_preco(false)
    { }

    PetscInt max_it; ///< Maximum iteration

    PetscReal rtol;  ///< Tolerance for the relative error, \f$e=|r|/|b|\f$.
    PetscReal atol;  ///< Tolerance for the absolute error, \f$e=|r|\f$.
    PetscReal dtol;  ///< Relative increase in the residual.

    /*!
        The name of solver, and it could be one of the following names
           "richardson"
           "chebychev"
           "cg"
           "cgne"
           "nash"
           "stcg"
           "gltr"
           "gmres"
           "fgmres"
           "lgmres"
           "dgmres"
           "tcqmr"
           "bcgs"
           "ibcgs"
           "bcgsl"
           "cgs"
           "tfqmr"
           "cr"
           "lsqr"
           "preonly"
           "qcg"
           "bicg"
           "minres"
           "symmlq"
           "lcd"
           "python"
           "broyden"
           "gcr"
           "ngmres"
           "specest"
     */
    std::string solver_name;

    /*!
         The name of preconditioner, and it could be one of the following names
          none
            jacobi
            sor
            lu
            shell
            bjacobi
            mg
            eisenstat
            ilu
            icc
            asm
            gasm
            ksp
            composite
            redundant
            spai
            nn
            cholesky
            pbjacobi
            mat
            hypre
            parms
            fieldsplit
            tfs
            ml
            prometheus
            galerkin
            exotic
            hmpi
            supportgraph
            asa
            cp
            bfbt
            lsc
            python
            pfmg
            syspfmg
            redistribute
            sacusp
            sacusppoly
            bicgstabcusp
            svd
            ainvcusp
            gamg
      */
    std::string _pc_name;

    /// Damping factor for Richards.
    double damping_factor_richards;

    /// Smallest eignvalue for Chebyshev.
    double emin_chebyshev;
    /// maximum eignvalue for Chebyshev.
    double emax_chebyshev;

    ///  Restart number of GMRES.
    double restart_number_gmres;

    /// Flag for Gram-Schmidt orthogonalization.
    bool is_gram_schmidt_orthogo_gmres

    /*!
         \brief Refinement type for GMRES.
         This iterative refinement is used to improve the stability
         of orthogonalization.
         KSPGMRESCGSRefinementType is a enum type of PETSc defined as
         typedef enum {KSP_GMRES_CGS_REFINE_NEVER, KSP_GMRES_CGS_REFINE_IFNEEDED,
                       KSP_GMRES_CGS_REFINE_ALWAYS} KSPGMRESCGSRefinementType;
    */
    KSPGMRESCGSRefinementType refine_type_gmres;

    /// Flag for right side preconditioning.
    bool righ_side_preco.
};

} // end namespace
#endif

