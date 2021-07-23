/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::invesionAlgorithms::PD

Description
    

SourceFiles
    PD.C

\*---------------------------------------------------------------------------*/

#ifndef PD_H
#define PD_H

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

class PD
:
    public inversionAlgorithm
{
    // Private data

public:

    //- Runtime type information
    TypeName("PD");


    // Constructors

        //- Construct from populationBalance and QuadratureMethod
        PD
        (
            const populationBalance& pb,
            PtrList<volScalarField>& nodes,
            PtrList<volScalarField>& weights
        );


    //- Destructor
    virtual ~PD();


    // Member Functions

        // Member Functions

        //- Calculates nodes and weights of the reconstructed NDF
        virtual void inversion(const Field<scalar>&);

        //- Calculates nodes and weights of the reconstructed NDF on the
        //- fixed value boundaries
        virtual void inversionBoundary();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace inversionAlgorithm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
