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

#include "qmom.H"
#include "inversionAlgorithm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace quadratureMethods
{
    defineTypeNameAndDebug(qmom, 0);
    addToRunTimeSelectionTable(quadratureMethod, qmom, populationBalance);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadratureMethods::qmom::qmom(const populationBalance& pb)
:
    quadratureMethod(pb)
{
    if (pb.numOfCoord() != 1)
    {
        FatalErrorInFunction
            << "QMOM is a valid option for univariate distributions:"
            << endl << "Number of coordinates has to be equal to 1"
            << endl << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::quadratureMethods::qmom::~qmom()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::quadratureMethods::qmom::nodesAndWeights
(
    const Field<scalar>& validCells
)
{
    // Using inversion algorithm to calculate nodes and weights
    invAlgm_->inversion(validCells);
    // Info<< "Nodes and weights calculated\n" << endl;
}


// ************************************************************************* //
