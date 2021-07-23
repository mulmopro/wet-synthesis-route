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
    Foam::BreakageKernels::powerLaw

Description
    

SourceFiles
    powerLawBreakage.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawBreakage_H
#define powerLawBreakage_H

#include "breakageKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakageKernels
{

/*---------------------------------------------------------------------------*\
                           Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    private breakageKernel
{
    // Private data

        //- Coefficient
        const scalar C1_;

        //- shear rate exponent
        const scalar C2_;

        //- diameter exponent
        const scalar C3_;

public:

    // Constructors

        //- Construct from a dictionary and a turbulence model
        powerLaw
        (
            const dictionary& dict,
            const incompressible::momentumTransportModel& turbulence
        );


    //- Destructor
    virtual ~powerLaw();


    //- Member Functions

        //- Breakage Frequency
        virtual tmp<DimensionedField<scalar, volMesh>> frequency
        (
            const DimensionedField<scalar, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakageKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "powerLawBreakage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
