/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 * \date   2012-05-10, 2014.10.10
 * \brief  Definition of the CPUTime class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <ctime>

namespace BaseLib
{
/// Count CPU time
class CPUTime
{
public:
    /// Start the timer.
    void start() { start_time_ = clock(); }

    /// Get the elapsed time after started.
    double elapsed() const
    {
        return (clock() - start_time_) / static_cast<double>(CLOCKS_PER_SEC);
    }

private:
    double start_time_ = 0.;
};

}  // end namespace BaseLib
