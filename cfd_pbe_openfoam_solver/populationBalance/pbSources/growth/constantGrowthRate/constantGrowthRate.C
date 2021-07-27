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

#include "constantGrowthRate.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace growthModels
{
    defineTypeNameAndDebug(constantGrowthRate, 0);
    addToRunTimeSelectionTable(growthModel, constantGrowthRate, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthModels::constantGrowthRate::constantGrowthRate
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    growthModel(dict, pb),

    superSat_(pb.superSat()),

    G0_("constGrowthRate", dimLength/dimTime, dict.lookup("G0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::growthModels::constantGrowthRate::~constantGrowthRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthModels::constantGrowthRate::rate
(
    const DimensionedField<scalar, volMesh>& Li
) const
{   
    tmp<DimensionedField<scalar, volMesh>> tx
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "constGrowthRate",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            G0_
        )
    );

    return tx*pos(superSat_ - 1.0);
}


realtype Foam::growthModels::constantGrowthRate::rate(realtype superSat) const
{
    if (superSat > 1.0)
    {
        return G0_.value();
    }

    return 0.0;
}

// ************************************************************************* //
