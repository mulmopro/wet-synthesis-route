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
    Foam::quadratureMethods::qmom

Description
    

SourceFiles
    qmom.C

\*---------------------------------------------------------------------------*/

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
