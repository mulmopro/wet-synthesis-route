/* ----------------------------------------------------------------------------
Copyright © 2021 Politecnico di Torino

This file is part of WetSynthRoute.

WetSynthRoute is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WetSynthRoute is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WetSynthRoute.  If not, see <https://www.gnu.org/licenses/>.

This offering is not approved or endorsed by OpenCFD Limited (ESI Group)
and the OpenFOAM Foundation.
---------------------------------------------------------------------------- */

#ifndef UserData_H
#define UserData_H

#include "sundials_types.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionNMC;
class populationBalance;
class growthModel;
class nucleationRateModel;
class nucleateSizeModel;

/*---------------------------------------------------------------------------*\
                      Struct userData Declaration
\*---------------------------------------------------------------------------*/

struct UserData{
    UserData
    (
        solutionNMC& solution,
        const growthModel& growth,
        const nucleationRateModel& nucRate,
        const nucleateSizeModel& nucleateSize
    ): solution_(solution), growth_(growth), nucRate_(nucRate),
    nucleateSize_(nucleateSize)
    {}
    label cell_id, nMetals, nMoments;
    realtype superSat, totalNH3, inertCharges, concNH3, concOH;
    realtype kv, crystalRho, crystalMW, effectiveConc;
    realtype metalCationConcs[3];
    solutionNMC& solution_;
    const growthModel& growth_;
    const nucleationRateModel& nucRate_;
    const nucleateSizeModel& nucleateSize_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
