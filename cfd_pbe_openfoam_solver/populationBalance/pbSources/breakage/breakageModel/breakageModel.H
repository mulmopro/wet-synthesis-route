/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::breakageModel_H

Description
    class for breakage phenomenon

SourceFiles
    breakageModel_H.C

\*---------------------------------------------------------------------------*/

#ifndef breakageModel_H
#define breakageModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "breakage.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class breakageModel Declaration
\*---------------------------------------------------------------------------*/

template<class BKernel, class BDaughterDist>
class breakageModel
:
    public breakage,
    public BKernel, 
    public BDaughterDist
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        breakageModel(const breakageModel&);

        //- Disallow default bitwise assignment
        void operator=(const breakageModel&);


public:

    //- Runtime type information
    TypeName("breakageModel");


    // Constructors

        //- Construct from components
        breakageModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~breakageModel()
    {}


    // Member Functions

        //- Return breakage source for moment k
        virtual tmp<Foam::DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "breakageModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
