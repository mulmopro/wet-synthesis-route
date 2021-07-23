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
    Foam::growthModels::constantGrowthRate

Description

SourceFiles
    constantGrowthRate.C

\*---------------------------------------------------------------------------*/

#ifndef constantGrowthRate_H
#define constantGrowthRate_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class populationBalance;

namespace growthModels
{

/*---------------------------------------------------------------------------*\
                           Class constantGrowthRate Declaration
\*---------------------------------------------------------------------------*/

class constantGrowthRate
:
    public growthModel
{
    //- Private members

    const DimensionedField<scalar, volMesh>& superSat_;
    
    dimensionedScalar G0_;

    
public:

    //- Runtime type information
    TypeName("constantRate");


    // Constructors

        //- Construct from a dictionary and a population balance model
        constantGrowthRate
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~constantGrowthRate();


    // Member Functions

        //- Constant growth rate
        tmp<DimensionedField<scalar, volMesh>>
            rate(const DimensionedField<scalar, volMesh>&) const;

        realtype rate(realtype superSat) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
