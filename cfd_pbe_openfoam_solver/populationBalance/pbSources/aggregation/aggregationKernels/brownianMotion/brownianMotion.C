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

#include "brownianMotion.H"
#include "kinematicMomentumTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aggregationKernels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

brownianMotion::brownianMotion
(
    const dictionary& dict,
    const incompressible::momentumTransportModel& turbulence
)
:
    aggregationKernel(turbulence),
    T_("T", dimTemperature, dict),
    rhoLiq_("rhoLiq", dimDensity, dict),
    CB_("CB", dimless, dict.lookupOrDefault<scalar>("CB", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

brownianMotion::~brownianMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
brownianMotion::frequency
(
    const DimensionedField<scalar, volMesh>& Li,
    const DimensionedField<scalar, volMesh>& Lj
) const
{
    return
        CB_
      * 2 * KB_ * T_ / 3.0
      / (turbulence_.nu()().internalField() * rhoLiq_)
      * (
            Li / max(Lj, dimensionedScalar("small", dimLength, SMALL))
          + Lj / max(Li, dimensionedScalar("small", dimLength, SMALL))
          + 2.0
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace Foam

// ************************************************************************* //
