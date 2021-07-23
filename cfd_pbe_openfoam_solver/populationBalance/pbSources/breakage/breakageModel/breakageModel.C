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

\*---------------------------------------------------------------------------*/

#include "breakageModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BKernel, class BDaughterDist>
Foam::breakageModel<BKernel, BDaughterDist>::breakageModel
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    breakage(dict, pb),
    BKernel(dict, pb.turbulence()),
    BDaughterDist(dict.subDict("daughterDistribution"), pb.turbulence())
    
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BKernel, class BDaughterDist>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::breakageModel<BKernel, BDaughterDist>::source
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
                "BSource" + Foam::name(k),
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

    forAll(nodes, i)
    {
        const DimensionedField<scalar, volMesh>& nodei
            = nodes[i].internalField();

        x +=
            BKernel::frequency(nodei)
          * weights[i].internalField()
          *
          (
              BDaughterDist::distribution(nodei, k)
            - Foam::pow(nodei, k)
          );
    }

    return tx;
}

// ************************************************************************* //
