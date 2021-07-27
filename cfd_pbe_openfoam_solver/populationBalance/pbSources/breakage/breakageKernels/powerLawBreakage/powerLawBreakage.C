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
    C1_(readScalar(dict.lookup("C1"))),
    C2_(readScalar(dict.lookup("C2"))),
    C3_(readScalar(dict.lookup("C3")))
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
    const fvMesh& mesh(Li.mesh());

    tmp<DimensionedField<scalar, volMesh>> tx
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "powerLawBreakageRate",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless/dimTime,
            C1_
          * pow
            (
                turbulence_.epsilon()().field() / turbulence_.nu()().field(),
                C2_ / 2.0
            )
          * pow(Li.field(), C3_)
        )
    );

    return tx;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakageKernels
} // End namespace Foam

// ************************************************************************* //
