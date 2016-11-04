/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   CapillaryPressureSaturation.h
 *
 */

#ifndef OGS_CAPILLARY_PRESSURE_SATURATION_H
#define OGS_CAPILLARY_PRESSURE_SATURATION_H

#include <string>

namespace MaterialLib
{
namespace PorousMedium
{
/// Base class of capillary pressure models
class CapillaryPressureSaturation
{
public:
    virtual ~CapillaryPressureSaturation() = default;

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get capillary pressure.
    virtual double getCapillaryPressure(const double saturation) const = 0;

    /// Get capillary pressure.
    virtual double getSaturation(const double capillary_ressure) const = 0;

    /// Get the derivative of the capillary pressure with respect to saturation
    virtual double getdPcdS(const double saturation) const = 0;

protected:
    /** A small number for an offset:
     *  1. to set the bound of S, the saturation, such that
     *     S in  [_Sr+_minor_offset, _Smax-_minor_offset]
     *  2. to set the bound of Pc, the capillary pressure, such that
     *     Pc in [_minor_offset, _Pc_max]
     */
    const double _minor_offset = std::numeric_limits<double>::epsilon();
};

}  // end namespace
}  // end namespace

#endif /* OGS_CAPILLARY_PRESSURE_SATURATION_H */
