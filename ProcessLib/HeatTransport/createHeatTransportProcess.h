/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createHeatTransportProcess.h
 *
 */

#ifndef OGS_CREATE_HEAT_TRANSPORT_PROCESS_H
#define OGS_CREATE_HEAT_TRANSPORT_PROCESS_H

#include <memory>
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace HeatTransport
{
std::unique_ptr<Process> createHeatTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);
}  // end of namespace
}  // end of namespace

#endif /* OGS_CREATE_HEAT_TRANSPORT_PROCESS_H */
