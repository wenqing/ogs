/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-12-06
 * \brief  Definition of the TemplateLogogFormatterSuppressedGCC class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATELOGOGFORMATTERSUPPRESSEDGCC_H_
#define TEMPLATELOGOGFORMATTERSUPPRESSEDGCC_H_

#include <string>

#if defined(USE_PETSC)
#include <petsc.h>
#elif defined(USE_MPI)
#include <mpi.h>
#endif

#include "logog/include/logog.hpp"

#include "StringTools.h"

namespace BaseLib {

/**
 * \brief TemplateLogogFormatterSuppressedGCC strips topics given as a template
 * parameter from logog::FormatterGCC.
 * See http://johnwbyrd.github.com/logog/customformatting.html for details.
 **/
template <int T_SUPPPRESS_TOPIC_FLAG>
class TemplateLogogFormatterSuppressedGCC : public logog::FormatterGCC
{
public:
#if defined(USE_MPI)
	TemplateLogogFormatterSuppressedGCC(MPI_Comm mpi_comm = MPI_COMM_WORLD);
#elif defined(USE_PETSC)
	TemplateLogogFormatterSuppressedGCC(MPI_Comm mpi_comm = PETSC_COMM_WORLD);
#endif

	virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
	{
	return ( logog::Formatter::GetTopicFlags( topic ) &
		~( T_SUPPPRESS_TOPIC_FLAG ));
	}

	virtual LOGOG_STRING &Format( const logog::Topic &topic, const logog::Target &target );

private:
#if defined(USE_MPI) || defined(USE_PETSC)
	std::string _str_mpi_rank;
#endif
};

} // namespace BaseLib

#include "TemplateLogogFormatterSuppressedGCC-impl.h"

#endif // TEMPLATELOGOGFORMATTERSUPPRESSEDGCC_H_
