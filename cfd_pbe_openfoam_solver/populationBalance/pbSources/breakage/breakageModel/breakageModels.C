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

\*---------------------------------------------------------------------------*/

#include "makeBreakageModels.H"

#include "constantRate.H"
#include "powerLawBreakage.H"

#include "erosion.H"
#include "Laakkonen.H"
#include "parabolic.H"
#include "symmetric.H"
#include "uniform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBreakageModel(constantRate, erosion);
makeBreakageModel(constantRate, Laakkonen);
makeBreakageModel(constantRate, parabolic);
makeBreakageModel(constantRate, symmetric);
makeBreakageModel(constantRate, uniform);

makeBreakageModel(powerLaw, erosion);
makeBreakageModel(powerLaw, Laakkonen);
makeBreakageModel(powerLaw, parabolic);
makeBreakageModel(powerLaw, symmetric);
makeBreakageModel(powerLaw, uniform);


