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

#include "constantGrowthRate.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace growthModels
{
    defineTypeNameAndDebug(constantGrowthRate, 0);
    addToRunTimeSelectionTable(growthModel, constantGrowthRate, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthModels::constantGrowthRate::constantGrowthRate
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    growthModel(dict, pb),

    superSat_(pb.superSat()),

    G0_("constGrowthRate", dimLength/dimTime, dict.lookup("G0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::growthModels::constantGrowthRate::~constantGrowthRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthModels::constantGrowthRate::rate
(
    const DimensionedField<scalar, volMesh>& Li
) const
{   
    tmp<DimensionedField<scalar, volMesh>> tx
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "constGrowthRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            G0_
        )
    );

    return tx*pos(superSat_ - 1.0);
}


realtype Foam::growthModels::constantGrowthRate::rate(realtype superSat) const
{
    if (superSat > 1.0)
    {
        return G0_.value();
    }

    return 0.0;
}

// ************************************************************************* //
