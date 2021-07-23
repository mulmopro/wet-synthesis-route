/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Struct
    Foam::userData

Description


\*---------------------------------------------------------------------------*/

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
    realtype kv, crystalRho, crystalMW;
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
