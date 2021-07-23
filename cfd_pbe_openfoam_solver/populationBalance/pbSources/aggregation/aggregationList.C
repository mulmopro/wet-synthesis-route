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

\*---------------------------------------------------------------------------*/

#include "aggregationList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aggregationList::aggregationList
(
        const dictionary& dict,
        const populationBalance& pb
)
:
    PtrList<aggregation>(),
    pb_(pb)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            Info<< "Creating '" << name << "' aggregation model:" << endl;

            this->set
            (
                i++,
                aggregation::New(modelDict, pb)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aggregationList::~aggregationList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::aggregationList::source
(
    const PtrList<volScalarField>& nodes,
    const PtrList<volScalarField>& weights,
    int k
) const
{
    tmp<DimensionedField<scalar, volMesh>> tx
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "CSource" + Foam::name(k),
                pb_.mesh().time().timeName(),
                pb_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pb_.mesh(),
            dimensionedScalar("zero", Foam::pow(dimLength, k - 3)/dimTime, 0)
        )
    );

    // Reference to tx in order to prevent referencing in the loop
    DimensionedField<scalar, volMesh>& x = tx.ref();

    forAll(*this, i)
    {
        x += operator[](i).source(nodes, weights, k);
    }

    return tx;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //



// ************************************************************************* //
