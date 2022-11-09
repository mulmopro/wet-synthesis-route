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

#include "erosion.H"
#include "kinematicMomentumTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace daughterDistributions
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

erosion::erosion
(
    const dictionary& dict,
    const incompressible::momentumTransportModel& turbulence
)
:
    daughterDistribution(turbulence),
    m_(readScalar(dict.lookup("m")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

erosion::~erosion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<DimensionedField<scalar, volMesh>>
erosion::distribution
(
    const DimensionedField<scalar, volMesh>& Li,
    int k
) const
{
    return
        pow(Li, k) * ((1 + pow(m_ - 1, k/3))/pow(m_, k/3));
}


scalar erosion::distribution(scalar Li, int k) const
{
    return
        pow(Li, k) * ((1 + pow(m_ - 1, k/3))/pow(m_, k/3));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterDistributions
} // End namespace Foam

// ************************************************************************* //
