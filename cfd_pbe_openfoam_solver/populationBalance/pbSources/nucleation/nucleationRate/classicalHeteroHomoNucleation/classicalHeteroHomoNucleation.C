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

#include "classicalHeteroHomoNucleation.H"
#include "populationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationRateModels
{
    defineTypeNameAndDebug(classicalHeteroHomoNucleation, 0);
    addToRunTimeSelectionTable(nucleationRateModel, classicalHeteroHomoNucleation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationRateModels::classicalHeteroHomoNucleation::classicalHeteroHomoNucleation
(
    const dictionary& dict,
    const populationBalance& pb
)
:
    nucleationRateModel(dict, pb),

    superSat_(pb.superSat()),
    
    rateDimension_("nucRateDim", dimless/dimVol/dimTime, 1.0),

    KJ_1_(readScalar(dict.lookup("KJ_1"))),

    KJ_2_(readScalar(dict.lookup("KJ_2"))),

    BJ_1_(readScalar(dict.lookup("BJ_1"))),

    BJ_2_(readScalar(dict.lookup("BJ_2")))
    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationRateModels::classicalHeteroHomoNucleation::~classicalHeteroHomoNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationRateModels::classicalHeteroHomoNucleation::rate()
const
{
  
    DimensionedField<scalar, volMesh> sqrLogS
    (
        Foam::sqr
        (
            Foam::log(max(superSat_, dimensionedScalar("minS", dimless, SMALL)))
        )
    );

    tmp<DimensionedField<scalar, volMesh>> tx
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "classicalHeteroHomoNucleation",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            rateDimension_
        )
    );

  return
    tx
  * (
        Foam::pow(10, KJ_1_) * exp(-1.0 * BJ_1_ / sqrLogS)
      + Foam::pow(10, KJ_2_) * exp(-1.0 * BJ_2_ / sqrLogS)
    )
  * pos(superSat_ - 1.0);
        
}


realtype Foam::nucleationRateModels::classicalHeteroHomoNucleation::rate
(realtype superSat) const
{
    if (superSat > 1.0)
    {
        realtype sqrLogS = Foam::sqr(Foam::log(superSat));


        return
            Foam::pow(10, KJ_1_) * exp(-1.0 * BJ_1_ / sqrLogS)
          + Foam::pow(10, KJ_2_) * exp(-1.0 * BJ_2_ / sqrLogS);
    }

    return 0.0;
}


// ************************************************************************* //
