/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "uniform.H"
#include "kinematicMomentumTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace daughterDistributions
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniform::uniform
(
    const dictionary& dict,
    const incompressible::momentumTransportModel& turbulence
)
:
    daughterDistribution(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

uniform::~uniform()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
uniform::distribution
(
    const DimensionedField<scalar, volMesh>& Li,
    int k
) const
{
    
    return
        (6.0 / (k + 3.0)) * pow(Li, k);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterDistributions
} // End namespace Foam

// ************************************************************************* //
