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
    Foam::aggregationKernels::turbulentCollision

Description
    

SourceFiles
    turbulentCollision.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentCollision_H
#define turbulentCollision_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                           Class turbulentCollision Declaration
\*---------------------------------------------------------------------------*/

class turbulentCollision
:
    private aggregationKernel
{
    // Private data

        //- Coefficient
        const dimensionedScalar CT_;


public:

    // Constructors

        //- Construct from a dictionary and a turbulence model
        turbulentCollision
        (
            const dictionary& dict,
            const incompressible::momentumTransportModel& turbulence
        );


    //- Destructor
    virtual ~turbulentCollision();


    //- Member Functions

        //- Aggregation Frequency
        virtual tmp<DimensionedField<scalar, volMesh>> frequency
        (
            const DimensionedField<scalar, volMesh>&,
            const DimensionedField<scalar, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
