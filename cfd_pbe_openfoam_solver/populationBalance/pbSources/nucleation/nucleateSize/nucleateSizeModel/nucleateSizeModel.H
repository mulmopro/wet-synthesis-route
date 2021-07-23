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
    Foam::nucleateSizeModel

Description

SourceFiles
    nucleateSizeModel.C
    newnucleateSizeModel.C

\*---------------------------------------------------------------------------*/

#ifndef nucleateSizeModel_H
#define nucleateSizeModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "sundials_types.h"

namespace Foam
{

class populationBalance;

/*---------------------------------------------------------------------------*\
                           Class nucleateSizeModel Declaration
\*---------------------------------------------------------------------------*/

class nucleateSizeModel
{
protected:

    // Protected data

    const fvMesh& mesh_; 

    const populationBalance& pb_;


public:

    //- Runtime type information
    TypeName("nucleateSizeModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleateSizeModel,
        dictionary,
        (
            const dictionary& dict,
            const populationBalance& pb
        ),
        (dict, pb)
    );


    // Constructors

        //- Construct from a dictionary and a population balance
        nucleateSizeModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~nucleateSizeModel();


    // Selectors

        static autoPtr<nucleateSizeModel> New
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    // Member Functions

        //- Nucleate size
        virtual tmp<DimensionedField<scalar, volMesh>> size() const = 0;

        virtual realtype size(realtype superSat) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
