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
    A_v_(A_.value()),
    rhoLiq_("rhoLiq", dimDensity, dict),
    rhoLiq_v_(rhoLiq_.value())
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


scalar bridgeFormation::efficiency
(
    scalar Li, scalar Lj, const PhysChemData& data
) const
{
    const scalar growthRate = data.growthRate;

    if (growthRate > SMALL)
    {
        const scalar epsilon = data.epsilon;

        const scalar nu = data.nu;

        scalar Li_pos = max(Li, SMALL);
        scalar Lj_pos = max(Lj, SMALL);

        scalar L_eq = Li*Lj / sqrt(pow(Li_pos - Lj_pos, 2) + Li_pos*Lj_pos);

        scalar Db =
            sqrt(rhoLiq_v_ / A_v_) * pow(epsilon*nu, 0.25) * L_eq;

        scalar r_L = max(Li_pos, Lj_pos) / min(Li_pos, Lj_pos);

        scalar sqrt_r_L = sqrt(r_L*r_L - 1.0);

        scalar f_lambda =
            4 * (1 + r_L - sqrt_r_L)
          / (
                (1.0/3.0 + r_L - sqrt_r_L)
              - pow(r_L - sqrt_r_L, 2) * (2*r_L / 3.0 + sqrt_r_L / 3.0)
            );

        return
            exp
            (
                -1.0
              * sqrt(epsilon / nu) * Db / f_lambda
              / max(growthRate, SMALL)
            );
    }

    return 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationEfficiencies
} // End namespace Foam

// ************************************************************************* //
