/*!
   \file
   \brief Definition of a function related to PETSc solver interface to assign
         the Dirichlet boundary conditions.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScTools.h"

namespace MathLib
{
void applyKnownSolution(PETScMatrix& A, PETScVector& b, PETScVector& x,
                        const std::vector<PetscInt>& vec_knownX_id,
                        const std::vector<PetscScalar>& vec_knownX_x)
{
    A.finalizeAssembly();

#ifndef NDEBUG
    if (std::any_of(vec_knownX_id.begin(), vec_knownX_id.end(),
                    [](PetscInt const i) { return i < 0; }))
    {
        OGS_FATAL(
            "Found negative indices in the vector of Dirichlet boundary "
            "conditions.");
    }
#endif  // NDEBUG

    A.setRowsColumnsZero(vec_knownX_id);
    A.finalizeAssembly();

    x.finalizeAssembly();
    b.finalizeAssembly();
    if (vec_knownX_id.size() > 0)
    {
        x.set(vec_knownX_id, vec_knownX_x);
        b.set(vec_knownX_id, vec_knownX_x);
    }

    x.finalizeAssembly();
    b.finalizeAssembly();
}

}  // end of namespace MathLib
