/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::breakage

Description

SourceFiles
    breakage.C

\*---------------------------------------------------------------------------*/

#ifndef breakage_H
#define breakage_H

#include "populationBalance.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class breakage Declaration
\*---------------------------------------------------------------------------*/

class breakage
{
protected:
    // Protected data

        //- Reference to the population balance
        const populationBalance& pb_;


public:

    //- Runtime type information
    TypeName("breakage");
    
    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            breakage,
            dictionary,
            (
                const dictionary& dict,
                const populationBalance& pb
            ),
            (dict, pb)
        );

    // Constructors

        //- Construct from the model table, dictionary and pairs
        breakage
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
        virtual ~breakage();

    // Selectors

        static autoPtr<breakage> New
        (
            const dictionary& dict,
            const populationBalance& pb
        );

    // Member Functions

        //- Return breakage source for moment k
        virtual tmp<DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
