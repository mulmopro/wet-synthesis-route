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
    Foam::growthModel

Description

SourceFiles
    growthModel.C
    newGrowthModel.C

\*---------------------------------------------------------------------------*/

#ifndef growthModel_H
#define growthModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "sundials_types.h"

namespace Foam
{

class populationBalance;

/*---------------------------------------------------------------------------*\
                           Class growthModel Declaration
\*---------------------------------------------------------------------------*/

class growthModel
{
protected:

    // Protected data

    const fvMesh& mesh_; 

    const populationBalance& pb_;


public:

    //- Runtime type information
    TypeName("growthModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        growthModel,
        dictionary,
        (
            const dictionary& dict,
            const populationBalance& pb
        ),
        (dict, pb)
    );


    // Constructors

        //- Construct from a dictionary and a population balance
        growthModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~growthModel();


    // Selectors

        static autoPtr<growthModel> New
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    // Member Functions

        //- Growth rate for each quadrature node
        virtual tmp<DimensionedField<scalar, volMesh>>
            rate(const DimensionedField<scalar, volMesh>&) const = 0;

        virtual realtype rate(realtype superSat) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
