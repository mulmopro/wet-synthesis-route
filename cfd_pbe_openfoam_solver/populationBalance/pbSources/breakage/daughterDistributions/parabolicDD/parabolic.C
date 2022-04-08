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

#include "parabolic.H"
#include "kinematicMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace daughterDistributions
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolic::parabolic
(
    const dictionary& dict,
    const incompressible::momentumTransportModel& turbulence
)
:
    daughterDistribution(turbulence),
    C_(readScalar(dict.lookup("C")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

parabolic::~parabolic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
parabolic::distribution
(
    const DimensionedField<scalar, volMesh>& Li,
    int k
) const
{
    
    return
        pow(Li, k)
      * (
            3 * C_ / (k + 3.0)
          + (1.0 - C_/2.0) * 18 * (6.0 - k) / ((k + 9.0)*(k + 6.0)*(k + 3.0))
        );
}


scalar parabolic::distribution(scalar Li, int k) const
{
    return
        pow(Li, k)
      * (
            3 * C_ / (k + 3.0)
          + (1.0 - C_/2.0) * 18 * (6.0 - k) / ((k + 9.0)*(k + 6.0)*(k + 3.0))
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterDistributions
} // End namespace Foam

// ************************************************************************* //
