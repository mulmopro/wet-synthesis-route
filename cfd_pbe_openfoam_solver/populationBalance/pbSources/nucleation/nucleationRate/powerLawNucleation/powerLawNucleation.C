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

#include "powerLawNucleation.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationRateModels
{
    defineTypeNameAndDebug(powerLawNucleation, 0);
    addToRunTimeSelectionTable(nucleationRateModel, powerLawNucleation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationRateModels::powerLawNucleation::powerLawNucleation
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    nucleationRateModel(dict, pb),

    superSat_(pb.superSat()),

    k_("nucRateCoeff", dimless/dimVol/dimTime, dict.lookup("k")),

    n_(readScalar(dict.lookup("n")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationRateModels::powerLawNucleation::~powerLawNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationRateModels::powerLawNucleation::rate()
const
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


realtype Foam::nucleationRateModels::powerLawNucleation::rate
(realtype superSat) const
{
    if (superSat > 1.0)
    {
        return k_.value() * Foam::pow((superSat - 1.0), n_);
    }

    return 0.0;
}


// ************************************************************************* //
