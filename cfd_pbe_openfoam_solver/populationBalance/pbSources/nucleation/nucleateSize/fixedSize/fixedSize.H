/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::nucleateSizeModels::fixedSize

Description

SourceFiles
    fixedSize.C

\*---------------------------------------------------------------------------*/

#ifndef fixedSize_H
#define fixedSize_H

#include "nucleateSizeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class populationBalance;

namespace nucleateSizeModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedSize Declaration
\*---------------------------------------------------------------------------*/

class fixedSize
:
    public nucleateSizeModel
{
    //- private members
    dimensionedScalar Xc_;


public:

    //- Runtime type information
    TypeName("fixed");


    // Constructors

        //- Construct from a dictionary and a population balance model
        fixedSize
        (
            const dictionary& dict,
            const populationBalance& pair
        );


    //- Destructor
    virtual ~fixedSize();


    // Member Functions

        //- Fixed nucleate size
        tmp<DimensionedField<scalar, volMesh>> size() const;

        realtype size(realtype superSat) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleateSizeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
