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

#include "activityCoeffModel.H"
#include "solutionNMC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(activityCoeffModel, 0);
    defineRunTimeSelectionTable(activityCoeffModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activityCoeffModel::activityCoeffModel
(
    const dictionary& dict,
    const solutionNMC& solution
)
:
    mesh_(solution.mesh()),
    solution_(solution)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activityCoeffModel::~activityCoeffModel()
{}


// ************************************************************************* //
