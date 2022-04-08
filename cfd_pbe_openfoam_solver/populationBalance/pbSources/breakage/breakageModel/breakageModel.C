/* ----------------------------------------------------------------------------
Copyright © 2021 Politecnico di Torino

This file is part of WetSynthRoute.

WetSynthRoute is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WetSynthRoute is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WetSynthRoute.  If not, see <https://www.gnu.org/licenses/>.

This offering is not approved or endorsed by OpenCFD Limited (ESI Group)
and the OpenFOAM Foundation.
---------------------------------------------------------------------------- */

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


template<class BKernel, class BDaughterDist> Foam::scalar
Foam::breakageModel<BKernel, BDaughterDist>::source
(
    const List<scalar>& nodes,
    const List<scalar>& weights,
    int k,
    const PhysChemData& data
) const
{
    scalar source_b = 0.0;

    forAll(nodes, i)
    {
        const scalar nodei = nodes[i];

        source_b +=
            BKernel::frequency(nodei, data)
          * weights[i]
          *
          (
              BDaughterDist::distribution(nodei, k)
            - Foam::pow(nodei, k)
          );
    }

    return source_b;
}

// ************************************************************************* //
