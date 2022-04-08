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

#include "turbulentCollision.H"
#include "kinematicMomentumTransportModel.H"
#include "UserData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aggregationKernels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentCollision::turbulentCollision
(
    const dictionary& dict,
    const incompressible::momentumTransportModel& turbulence
)
:
    aggregationKernel(turbulence),
    CT_("CT", dimless, dict),
    CT_v_(CT_.value())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

turbulentCollision::~turbulentCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
turbulentCollision::frequency
(
    const DimensionedField<scalar, volMesh>& Li,
    const DimensionedField<scalar, volMesh>& Lj
) const
{
    const DimensionedField<scalar, volMesh> epsilon =
        turbulence_.epsilon()().internalField();
    // dimensionedScalar epsilon =
    //     dimensionedScalar("constEps", dimVelocity*dimVelocity/dimTime, 8.5);

    return
        CT_
      * 2.2943
      * sqrt
        (
            epsilon
          / turbulence_.nu()().internalField()
        )
      * pow(Li + Lj, 3);
}

scalar turbulentCollision::frequency
(
    scalar Li, scalar Lj, const PhysChemData& data
) const
{

    return
        CT_v_
      * 2.2943
      * sqrt
        (
            data.epsilon
          / data.nu
        )
      * pow(Li + Lj, 3);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace Foam

// ************************************************************************* //
