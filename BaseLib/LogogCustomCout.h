/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LOGOGCUSTOMCOUT_H_
#define LOGOGCUSTOMCOUT_H_

#include <ostream>

#if defined(USE_PETSC)
#include <petsc.h>
#elif defined(USE_MPI)
#include <mpi.h>
#endif

#include "logog/include/logog.hpp"

namespace BaseLib
{

/// Custom target for logog output
class LogogCustomCout : public logog::Target
{
public:
#if defined(USE_MPI) || defined(USE_PETSC)
	/**
	 * Constructor when MPI is involved
	 *
	 * @param all_rank_output_level  Minimum level to output messages from all MPI processes
	 * @param mpi_comm               MPI communicator
	 */
#if defined(USE_PETSC)
	LogogCustomCout(LOGOG_LEVEL_TYPE all_rank_output_level = LOGOG_LEVEL_INFO, MPI_Comm mpi_comm = PETSC_COMM_WORLD)
#elif defined(USE_MPI)
	LogogCustomCout(LOGOG_LEVEL_TYPE all_rank_output_level = LOGOG_LEVEL_INFO, MPI_Comm mpi_comm = MPI_COMM_WORLD)
#endif
	: _all_rank_output_level(all_rank_output_level), _is_rank0 (getRank(mpi_comm)==0)
	{}
#endif

	virtual int Receive( const logog::Topic &topic )
	{
#if defined(USE_MPI) || defined(USE_PETSC)
		if (topic.Level() > _all_rank_output_level && !_is_rank0)
			return 0;
#endif
		return logog::Target::Receive(topic);
	}

	virtual int Output( const LOGOG_STRING &data )
	{
#if defined(USE_MPI) || defined(USE_PETSC)
		if (_is_rank0)
#endif
		LOGOG_COUT << (const LOGOG_CHAR *)data << std::flush;
		return 0;
	}

private:
#if defined(USE_MPI) || defined(USE_PETSC)
	int getRank(MPI_Comm mpi_comm) const
	{
		int rank = 0;
		MPI_Comm_rank(mpi_comm, &rank);
		return rank;
	}

	const LOGOG_LEVEL_TYPE _all_rank_output_level;
	const bool _is_rank0;
#endif
};

#endif // LOGOGCUSTOMCOUT_H_

} // namespace BaseLib
