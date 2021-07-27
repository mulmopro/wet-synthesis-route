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

#ifndef adaptiveWheeler_H
#define adaptiveWheeler_H

#include "inversionAlgorithm.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace inversionAlgorithms
{

/*---------------------------------------------------------------------------*\
                           Class PD Declaration
\*---------------------------------------------------------------------------*/

class adaptiveWheeler
:
    public inversionAlgorithm
{
    // Private data

public:

    //- Runtime type information
    TypeName("adaptiveWheeler");


    // Constructors

        //- Construct from populationBalance and QuadratureMethod
        adaptiveWheeler
        (
            const populationBalance& pb,
            PtrList<volScalarField>& nodes,
            PtrList<volScalarField>& weights
        );


    //- Destructor
    virtual ~adaptiveWheeler();


    // Member Functions

        // Member Functions

        //- Calculates nodes and weights of the reconstructed NDF
        virtual void inversion(const Field<scalar>&);

        //- Calculates nodes and weights of the reconstructed NDF on the
        //- fixed value boundaries
        virtual void inversionBoundary();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace inversionAlgorithms
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
