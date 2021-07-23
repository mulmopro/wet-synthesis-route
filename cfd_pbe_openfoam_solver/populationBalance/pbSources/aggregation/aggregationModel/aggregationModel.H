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
    Foam::aggregationModel

Description
    Class for the aggregation model including both the kernel and efficiency

SourceFiles
    aggregationModel.C

\*---------------------------------------------------------------------------*/

#ifndef aggregationModel_H
#define aggregationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aggregation.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class aggregationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CKernel, class CEfficiency>
class aggregationModel
:
    public aggregation,
    public CKernel, 
    public CEfficiency
{

public:

    //- Runtime type information
    TypeName("aggregationModel");


    // Constructors

        //- Construct from components
        aggregationModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
        virtual ~aggregationModel()
        {}


    // Member Functions

        //- Return aggregation source for moment k
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
    #include "aggregationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
