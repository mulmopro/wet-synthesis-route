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

#include "powerLawGrowth.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace growthModels
{
    defineTypeNameAndDebug(powerLawGrowth, 0);
    addToRunTimeSelectionTable(growthModel, powerLawGrowth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::growthModels::powerLawGrowth::powerLawGrowth
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    growthModel(dict, pb),

    superSat_(pb.superSat()),

    k_("GrowthBySize", dimLength/dimTime, dict.lookup("k")),

    n_(readScalar(dict.lookup("n")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::growthModels::powerLawGrowth::~powerLawGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::growthModels::powerLawGrowth::rate
(
    const DimensionedField<scalar, volMesh>& Li
) const
{
    
    return
        k_
      * Foam::pow
        (
            max
            (
                (superSat_ - 1.0),
                dimensionedScalar("zero", dimless, 0.0)
            ),
            n_
        );
}


realtype Foam::growthModels::powerLawGrowth::rate(realtype superSat) const
{
    if (superSat > 1.0)
    {
        return k_.value() * Foam::pow((superSat - 1.0), n_);
    }

    return 0.0;
}


// ************************************************************************* //
