/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::aggregationList

Description
    List container for aggregation models

SourceFiles
    aggregationList.C

\*---------------------------------------------------------------------------*/

#ifndef aggregationList_H
#define aggregationList_H

#include "aggregation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aggregationList Declaration
\*---------------------------------------------------------------------------*/

class aggregationList
:
    PtrList<aggregation>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        aggregationList(const aggregationList&);

        //- Disallow default bitwise assignment
        void operator=(const aggregationList&);


protected:

    // Protected data

        //- Reference to the population balance
        const populationBalance& pb_;


public:

    //- Constructor
    aggregationList
    (
        const dictionary& dict,
        const populationBalance& pb
    );

    //- Destructor
    ~aggregationList();


    // Member Functions

        //- Calculate the aggregation models
        tmp<DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
