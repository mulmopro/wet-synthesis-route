/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "powerLawGrowth.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace growthModels
{
    defineTypeNameAndDebug(powerLawGrowth, 0);
    addToRunTimeSelectionTable(growthModel, powerLawGrowth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthModels::powerLawGrowth::powerLawGrowth
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    growthModel(dict, pb),

    superSat_(pb.superSat()),

    k_("GrowthBySize", dimLength/dimTime, dict.lookup("k")),

    n_(readScalar(dict.lookup("n")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::growthModels::powerLawGrowth::~powerLawGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthModels::powerLawGrowth::rate
(
    const DimensionedField<scalar, volMesh>& Li
) const
{
    
    return
        k_
      * Foam::pow
        (
            max
            (
                (superSat_ - 1.0),
                dimensionedScalar("zero", dimless, 0.0)
            ),
            n_
        );
}


realtype Foam::growthModels::powerLawGrowth::rate(realtype superSat) const
{
    if (superSat > 1.0)
    {
        return k_.value() * Foam::pow((superSat - 1.0), n_);
    }

    return 0.0;
}


// ************************************************************************* //
