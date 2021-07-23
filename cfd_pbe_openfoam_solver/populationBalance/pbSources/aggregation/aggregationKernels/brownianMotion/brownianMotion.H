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
    Foam::aggregationKernels::brownianMotion

Description
    

SourceFiles
    brownianMotion.C

\*---------------------------------------------------------------------------*/

#ifndef brownianMotion_H
#define brownianMotion_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                           Class brownianMotion Declaration
\*---------------------------------------------------------------------------*/

class brownianMotion
:
    private aggregationKernel
{
    // Private data

        //- Coefficient
        const dimensionedScalar KB_ =
            dimensionedScalar("KB", dimEnergy/dimTemperature, 1.38064852e-23);

        const dimensionedScalar T_;

        const dimensionedScalar rhoLiq_;

        const dimensionedScalar CB_;


public:

    // Constructors

        //- Construct from a dictionary and a turbulence model
        brownianMotion
        (
            const dictionary& dict,
            const incompressible::momentumTransportModel& turbulence
        );


    //- Destructor
    virtual ~brownianMotion();


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
    #include "brownianMotion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
