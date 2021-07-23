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

#include "bridgeFormation.H"
#include "populationBalance.H"
#include "growthModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aggregationEfficiencies
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bridgeFormation::bridgeFormation
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    aggregationEfficiency(pb),
    turbulence_(pb.turbulence()),
    growth_(pb.growth()),
    A_("A", dimForce/dimArea, dict),
    rhoLiq_("rhoLiq", dimDensity, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bridgeFormation::~bridgeFormation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
bridgeFormation::efficiency
(
    const DimensionedField<scalar, volMesh>& Li,
    const DimensionedField<scalar, volMesh>& Lj
) const
{
    const DimensionedField<scalar, volMesh> epsilon =
        turbulence_.epsilon()().internalField();
    // dimensionedScalar epsilon =
    //     dimensionedScalar("constEps", dimVelocity*dimVelocity/dimTime, 8.5);

    const DimensionedField<scalar, volMesh> nu =
        turbulence_.nu()().internalField();

    DimensionedField<scalar, volMesh> Li_pos =
        max(Li, dimensionedScalar("small", dimLength, SMALL));
    DimensionedField<scalar, volMesh> Lj_pos =
        max(Lj, dimensionedScalar("small", dimLength, SMALL));

    DimensionedField<scalar, volMesh> L_eq =
        Li*Lj / sqrt(pow(Li_pos - Lj_pos, 2) + Li_pos*Lj_pos);

    DimensionedField<scalar, volMesh> growthRate = growth_.rate(L_eq);

    DimensionedField<scalar, volMesh> Db =
        sqrt(rhoLiq_ / A_) * pow(epsilon*nu, 0.25) * L_eq;

    DimensionedField<scalar, volMesh> r_L =
        max(Li_pos, Lj_pos) / min(Li_pos, Lj_pos);

    DimensionedField<scalar, volMesh> sqrt_r_L = sqrt(r_L*r_L - 1.0);

    DimensionedField<scalar, volMesh> f_lambda =
        4 * (1 + r_L - sqrt_r_L)
      / (
            (1.0/3.0 + r_L - sqrt_r_L)
          - pow(r_L - sqrt_r_L, 2) * (2*r_L / 3.0 + sqrt_r_L / 3.0)
        );

    return
        pos(growthRate)
      * exp
        (
            -1.0
          * sqrt(epsilon / nu) * Db / f_lambda
          / max(growthRate, dimensionedScalar("small", dimVelocity, SMALL))
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationEfficiencies
} // End namespace Foam

// ************************************************************************* //
