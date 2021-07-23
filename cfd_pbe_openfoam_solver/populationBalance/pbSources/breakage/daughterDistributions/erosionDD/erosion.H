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
    Foam::daughterDistributions::erosion

Description
    

SourceFiles
    erosion.C

\*---------------------------------------------------------------------------*/

#ifndef erosionDD_H
#define erosionDD_H

#include "daughterDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace daughterDistributions
{

/*---------------------------------------------------------------------------*\
                         Class erosion Declaration
\*---------------------------------------------------------------------------*/

class erosion
:
    public daughterDistribution
{
    // Private data

public:

    // Constructors

        //- Construct from a dictionary and a phase pair
        erosion
        (
            const dictionary& dict,
            const incompressible::momentumTransportModel& turbulence
        );


    //- Destructor
    virtual ~erosion();


    // Member Functions

        //- Daughter distribution divided by Li to power k
         virtual tmp<DimensionedField<scalar, volMesh>> distribution
        (
            const DimensionedField<scalar, volMesh>&,
            int
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterDistributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "erosion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
