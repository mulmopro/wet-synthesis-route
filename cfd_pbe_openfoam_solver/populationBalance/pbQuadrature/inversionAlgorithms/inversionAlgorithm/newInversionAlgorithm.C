/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "inversionAlgorithm.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::inversionAlgorithm> Foam::inversionAlgorithm::New
(
    const populationBalance& pb,
    PtrList<volScalarField>& nodes,
    PtrList<volScalarField>& weights
)
{
    word inversionAlgorithmType
    (
        pb.subDict("quadrature").lookup("inversionAlgorithm")
    );

    Info<< "Selecting inversion algorithm for "
        << "quadrature method" << ": " << inversionAlgorithmType << endl;

    populationBalanceConstructorTable::iterator cstrIter =
        populationBalanceConstructorTablePtr_->find(inversionAlgorithmType);

    if (cstrIter == populationBalanceConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown inverison algorithm "
            << inversionAlgorithmType << endl << endl
            << "Valid inversion algorithms are : " << endl
            << populationBalanceConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(pb, nodes, weights);
}


// ************************************************************************* //
