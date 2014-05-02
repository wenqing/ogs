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
    PETSc KSP Richards options

 */
struct PETScPC_KSP_Richards_Option
{
    PETScPC_KSP_Richards_Option(const boost::property_tree::ptree &option);

    /// Set Richards option
    void setOption(KSP &solver)
    {
        KSPRichardsonSetScale(solver, damping_factor_richards);
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_KSP_Richards_Option& opt)
    {
        damping_factor_richards = opt.damping_factor_richards;
    }

    /// Damping factor for Richards.
    PetscReal damping_factor_richards;
};

/*!
    PETSc KSP Chebyshev options

 */
struct PETScPC_KSP_Chebyshev_Option
{
    PETScPC_KSP_Chebyshev_Option(const boost::property_tree::ptree &option);

    /// Set Chebyshev option
    void setOption(KSP &solver)
    {
        KSPChebyshevSetEigenvalues(solver, emax_chebyshev, emin_chebyshev);
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_KSP_Chebyshev_Option& opt)
    {
        emin_chebyshev = opt.emin_chebyshev;
        emax_chebyshev = opt.emax_chebyshev;
    }

    /// Smallest eignvalue for Chebyshev.
    PetscReal emin_chebyshev;
    /// maximum eignvalue for Chebyshev.
    PetscReal emax_chebyshev;
};

/*!
    PETSc KSP GMRES options

 */
struct PETScPC_KSP_GMRES_Option
{
    PETScPC_KSP_GMRES_Option(const boost::property_tree::ptree &option);

    /// Set GMRES option
    void setOption(KSP &solver);

    /// Overloaded assign operator
    void operator = (const PETScPC_KSP_GMRES_Option& opt)
    {
        restart_number_gmres = opt.restart_number_gmres;
        is_modified_gram_schmidt_gmres = opt.is_modified_gram_schmidt_gmres;
        refine_type_gmres = opt.refine_type_gmres;
    }

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
};

/*!
    PETSc ILU preconditioner options

 */
struct PETScPC_ILU_Option
{
    PETScPC_ILU_Option(const boost::property_tree::ptree &option);

    /// Set ILU option
    void setOption(PC &pc);

    /// Overloaded assign operator
    void operator = (const PETScPC_ILU_Option& opt)
    {
        levels = opt.levels;
        reuse_ordering = opt.reuse_ordering;
        reuse_fill = opt.reuse_fill;
        use_in_place = opt.use_in_place;
        allow_diagonal_fill = opt.allow_diagonal_fill;
    }

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
    PETScPC_SOR_Option(const boost::property_tree::ptree &option);

    /// Set SOR/SSOR option
    void setOption(PC &pc);

    /// Overloaded assign operator
    void operator = (const PETScPC_SOR_Option& opt)
    {
        omega = opt.omega;
        its = opt.its;
        lits = opt.lits;
        type = opt.type;
    }

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
    PETSc LU preconditioner options

*/
struct PETScPC_LU_Option
{
    PETScPC_LU_Option(const boost::property_tree::ptree &option);

    /// Set LU option
    void setOption(PC &pc)
    {
        PCFactorSetMatOrderingType(pc, mat_type);
    }

    /// Overloaded assign operator
    void operator = (const PETScPC_LU_Option& opt)
    {
        mat_type = opt.mat_type;
    }

    /// Ordering routine (to reduce fill) to be used in the LU factorization.
    MatOrderingType mat_type;
};

/*!
   \brief Class to manage the basic configuration data for
          a Krylov subspace (KSP) interation method and a procondtioner (PC)
          of PETSc routines.
*/
class PETScLinearSolverOption
{
    public:

        /*!
            Default constructor
            \param ksp_option ptree of basic configuration data for solver
            \param pc_option  ptree of basic configuration data for preconditioner
        */
        PETScLinearSolverOption(const boost::property_tree::ptree &ksp_option,
                                const boost::property_tree::ptree &pc_option);

        /// Set basic options for a KSP and a PC
        void setOption(KSP &ksp, PC &pc)
        {
            KSPSetType(ksp,  _solver_name.c_str());
            KSPGetPC(ksp, &pc);
            PCSetType(pc, _pc_name.c_str());
            KSPSetPCSide(ksp, _preco_side);

            KSPSetTolerances(ksp, _rtol, _atol, _dtol, _max_it);
        }

    private:
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
        std::string _solver_name;

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

        /// Flag for which side preconditioning.
        PCSide _preco_side;

        PetscInt _max_it; ///< Maximum iteration.

        PetscReal _rtol;  ///< Tolerance for the relative error, \f$e=|r|/|b|\f$.
        PetscReal _atol;  ///< Tolerance for the absolute error, \f$e=|r|\f$.
        PetscReal _dtol;  ///< Relative increase in the residual.
};

} // end namespace
#endif

