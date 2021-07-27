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

#include "inversionAlgorithm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inversionAlgorithm, 0);
    defineRunTimeSelectionTable(inversionAlgorithm, populationBalance);
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * //

bool Foam::inversionAlgorithm::nodeReduction
(
    const Foam::List<scalar>& zeta,
    int& n
)
{
    bool realizableMoms = true;

    label counter = 1;

    if (zeta[0] <= 0)
    {
        n = 1;
        return false;
    }
    
    for (int i=2; i<2*n-1; i+=2)
    {
        if (zeta[i] <= 0 || zeta[i-1] <= 0)
        {
            n = counter;
            realizableMoms = false;
            break;
        }
        counter++;
    }

    return realizableMoms;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inversionAlgorithm::inversionAlgorithm
(
    const populationBalance& pb,
    PtrList<volScalarField>& nodes,
    PtrList<volScalarField>& weights
)
:
    moments_(pb.moments()),
    dSmall_(pb.dSmall()),
    numOfNodes_(pb.numOfNodes()),
    nodes_(nodes),
    weights_(weights)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inversionAlgorithm::~inversionAlgorithm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
