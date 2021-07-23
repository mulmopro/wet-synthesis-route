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

#include "aggregationModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CKernel, class CEfficiency>
Foam::aggregationModel<CKernel, CEfficiency>::aggregationModel
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    aggregation(dict, pb),
    CKernel(dict, pb.turbulence()),
    CEfficiency(dict.subDict("efficiency"), pb)
    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CKernel, class CEfficiency>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::aggregationModel<CKernel, CEfficiency>::source
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
                "CModelSource" + Foam::name(k),
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
        const DimensionedField<scalar, volMesh>& weighti
            = weights[i].internalField();

        const DimensionedField<scalar, volMesh>& nodei
            = nodes[i].internalField();

        const DimensionedField<scalar, volMesh> nodeiToPow3
        (
            Foam::pow(nodei, 3.0)
        );

        const DimensionedField<scalar, volMesh> nodeiToPowK
        (
            Foam::pow(nodei, k)
        );

        forAll(nodes, j)
        {
            const DimensionedField<scalar, volMesh>& nodej
                = nodes[j].internalField();

            x +=
                CKernel::frequency(nodei, nodej)
              * CEfficiency::efficiency(nodei, nodej)
              * weighti * weights[j].internalField()
              *
              (
                    0.5
                  * Foam::pow
                    (
                        nodeiToPow3
                      + Foam::pow(nodej, 3.0)
                      , k/3.0
                    )
                  - nodeiToPowK
              );
        }
    }

    return tx;
}

// ************************************************************************* //
