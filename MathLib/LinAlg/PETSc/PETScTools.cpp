/*!
   \file  PETScTools.cpp
   \brief Definition of a function related to PETSc solver interface to assign
         the Dirichlet boundary conditions.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScTools.h"

namespace MathLib
{

void applyKnownSolution(PETScMatrix &A, PETScVector &b, PETScVector &x,
                        const std::vector<PetscInt> &_vec_knownX_id,
                        const std::vector<PetscScalar> &vec_knownX_x);
{
    const PetscInt ni = static_cast<PetscInt> (_vec_knownX_id.size());

    A.setRowsColumnsZero(ni, &_vec_knownX_id[0]);
    A.finalizeAssembly;

    x.set<vector<PetscScalar>>(vec_knownX_id, vec_knownX_x);
    b.set<vector<PetscScalar>>(vec_knownX_id, vec_knownX_x);

    x.finalizeAssembly();
    b.finalizeAssembly();
}

} // end of namespace MathLib


