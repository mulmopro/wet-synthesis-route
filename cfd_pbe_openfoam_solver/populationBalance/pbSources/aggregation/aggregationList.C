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

#include "aggregationList.H"
#include "UserData.H"

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


Foam::scalar Foam::aggregationList::source
(
    const List<scalar>& nodes,
    const List<scalar>& weights,
    int k,
    const PhysChemData& data
) const
{
    scalar source_a = 0;

    forAll(*this, i)
    {
        source_a += operator[](i).source(nodes, weights, k, data);
    }

    return source_a;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //



// ************************************************************************* //
