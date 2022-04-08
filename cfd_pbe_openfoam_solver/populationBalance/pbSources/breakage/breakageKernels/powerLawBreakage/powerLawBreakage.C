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

#include "powerLawBreakage.H"
#include "kinematicMomentumTransportModel.H"
#include "UserData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakageKernels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

powerLaw::powerLaw
(
    const dictionary& dict,
    const incompressible::momentumTransportModel& turbulence
)
:
    breakageKernel(turbulence),
    Cbr_("Cbr", dimless, dict),
    Cbr_v_(Cbr_.value()),
    gamma_(readScalar(dict.lookup("gamma")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
powerLaw::frequency
(
    const DimensionedField<scalar, volMesh>& Li
) const
{
    const DimensionedField<scalar, volMesh> epsilon =
        turbulence_.epsilon()().internalField();

    const DimensionedField<scalar, volMesh> nu =
        turbulence_.nu()().internalField();

    return
        Cbr_
      * pow
        (
            Li / pow(pow(nu, 3) / epsilon, 0.25),
            gamma_
        )
      / pow(nu / epsilon, 0.5);
}


scalar powerLaw::frequency(scalar Li, const PhysChemData& data) const
{
    const scalar epsilon = data.epsilon;

    const scalar nu = data.nu;

    return
        Cbr_v_
      * pow
        (
            Li / pow(pow(nu, 3) / epsilon, 0.25),
            gamma_
        )
      / pow(nu / epsilon, 0.5);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakageKernels
} // End namespace Foam

// ************************************************************************* //
