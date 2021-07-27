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

#ifndef qmom_H
#define qmom_H

#include "quadratureMethod.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace quadratureMethods
{

/*---------------------------------------------------------------------------*\
                           Class qmom Declaration
\*---------------------------------------------------------------------------*/

class qmom
:
    public quadratureMethod
{
    // Private data

public:

    //- Runtime type information
    TypeName("QMOM");


    // Constructors

        //- Construct from a dictionary and a phase pair
        qmom(const populationBalance& pb);


    //- Destructor
    virtual ~qmom();


    // Member Functions

        //- Calculates nodes and weights of the reconstructed NDF
        virtual void nodesAndWeights(const Field<scalar>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace quadratureMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
