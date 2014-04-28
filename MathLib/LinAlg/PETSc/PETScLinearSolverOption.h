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
#include <boost/property_tree/ptree.hpp>

#include <petscksp.h>

namespace MathLib
{

/*!
    PETSc ILU preconditioner options

 */
struct PETScPC_ILU_Option
{
    PETScPC_ILU_Option();

    void setOption(const boost::property_tree::ptree &option);

    /// Number of levels of fill for ILU(k)
    int levels;

    /// Reuse ordering of factorized matrix from previous factorization
    bool reuse_ordering;

    /// Use the fill ratio computed in the initial factorization.
    bool reuse_fill;

    /*! for ILU(0) with natural ordering, reuses the space of the matrix
        for its factorization (overwrites original matrix)
    */
    bool use_in_place;

    /// fill in a zero diagonal even if levels of fill indicate it wouldn't be fill
    bool allow_diagonal_fill;
};

/*!
    PETSc SOR//SSOR preconditioner options

*/
struct PETScPC_SOR_Option
{
    PETScPC_SOR_Option();

    void setOption(const boost::property_tree::ptree &option);

    /// Relaxation number
    PetscReal omega;

    /// Number of parelllel iterations, each parallel iteration
    /// has 'lits' local iterations
    PetscInt its;

    /// Number of local iterations
    PetscInt lits;

    /*!
        SOR type:
        SOR_FORWARD_SWEEP
        SOR_BACKWARD_SWEEP
        SOR_SYMMETRIC_SWEEP
        SOR_LOCAL_FORWARD_SWEEP
        SOR_LOCAL_BACKWARD_SWEEP
        SOR_LOCAL_SYMMETRIC_SWEEP
    */
    MatSORType type;
};

/*!
   \brief This a struct data containing configuration data
          used to configure a procondtioner and a linear solver of PETSc.
*/
struct PETScLinearSolverOption
{
    /// Default constructor
    PETScLinearSolverOption()
        : solver_name("bcgs"), pc_name("bjacobi"), preco_side(PC_LEFT),
          max_it(2000), rtol(1.e-5), atol(PETSC_DEFAULT), dtol(PETSC_DEFAULT),
          damping_factor_richards(1.0), emin_chebyshev(0.01), emax_chebyshev(100.0),
          restart_number_gmres(30), is_modified_gram_schmidt_gmres(false),
          refine_type_gmres(KSP_GMRES_CGS_REFINE_NEVER)
    { }

    /// Constructor with an argument of options
    PETScLinearSolverOption(const boost::property_tree::ptree &option);

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
    std::string pc_name;

    /// Flag for which side preconditioning.
    PCSide preco_side;

    PetscInt max_it; ///< Maximum iteration

    PetscReal rtol;  ///< Tolerance for the relative error, \f$e=|r|/|b|\f$.
    PetscReal atol;  ///< Tolerance for the absolute error, \f$e=|r|\f$.
    PetscReal dtol;  ///< Relative increase in the residual.

    /// Damping factor for Richards.
    PetscReal damping_factor_richards;

    /// Smallest eignvalue for Chebyshev.
    PetscReal emin_chebyshev;
    /// maximum eignvalue for Chebyshev.
    PetscReal emax_chebyshev;

    ///  Restart number of GMRES.
    PetscInt restart_number_gmres;

    /// Flag for the modified Gram-Schmidt orthogonalization.
    bool is_modified_gram_schmidt_gmres;

    /*!
         \brief Refinement type for GMRES.
         This iterative refinement is used to improve the stability
         of orthogonalization.
         KSPGMRESCGSRefinementType is a enum type of PETSc defined as
         typedef enum {KSP_GMRES_CGS_REFINE_NEVER, KSP_GMRES_CGS_REFINE_IFNEEDED,
                       KSP_GMRES_CGS_REFINE_ALWAYS} KSPGMRESCGSRefinementType;
    */
    KSPGMRESCGSRefinementType refine_type_gmres;

    /// ilu or icc preconditioner options
    PETScPC_ILU_Option pc_ilu;

};

} // end namespace
#endif

