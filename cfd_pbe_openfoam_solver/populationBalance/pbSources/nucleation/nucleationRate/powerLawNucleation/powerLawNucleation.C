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

#include "powerLawNucleation.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationRateModels
{
    defineTypeNameAndDebug(powerLawNucleation, 0);
    addToRunTimeSelectionTable(nucleationRateModel, powerLawNucleation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationRateModels::powerLawNucleation::powerLawNucleation
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    nucleationRateModel(dict, pb),

    superSat_(pb.superSat()),

    k_("nucRateCoeff", dimless/dimVol/dimTime, dict.lookup("k")),

    n_(readScalar(dict.lookup("n")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationRateModels::powerLawNucleation::~powerLawNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationRateModels::powerLawNucleation::rate()
const
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


realtype Foam::nucleationRateModels::powerLawNucleation::rate
(realtype superSat) const
{
    if (superSat > 1.0)
    {
        return k_.value() * Foam::pow((superSat - 1.0), n_);
    }

    return 0.0;
}


// ************************************************************************* //
