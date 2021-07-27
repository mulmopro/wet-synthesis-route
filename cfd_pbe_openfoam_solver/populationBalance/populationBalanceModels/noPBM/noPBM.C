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

#include "noPBM.H"
#include "addToRunTimeSelectionTable.H"
#include "growthModel.H"
#include "nucleationRateModel.H"
#include "nucleateSizeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalances
{
    defineTypeNameAndDebug(noPBM, 0);
    addToRunTimeSelectionTable(populationBalance, noPBM, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalances::noPBM::noPBM
(
    const fvMesh& mesh,
    const incompressible::momentumTransportModel& turbulence,
    const solutionNMC& solution
)
:
    populationBalance(mesh, turbulence, solution),

    dSmall_(0),

    numOfNodes_(0),
    
    numOfMoments_(0),
    
    numOfCoord_(0),

    precRate_
    (
        IOobject
        (
            "precRate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMoles/dimVol/dimTime, 0)
    ),

    moments_(0),

    cationConcRatios_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalances::noPBM::~noPBM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::growthModel& Foam::populationBalances::noPBM::growth()
const
{
    FatalErrorInFunction
        << "Invalid request of growth model when no population balance model"
        << "is active" << exit(FatalError);

    return growthModel::New(*this, *this)();
}


const Foam::nucleationRateModel&
Foam::populationBalances::noPBM::nucleationRate()
const
{
    FatalErrorInFunction
        << "Invalid request of growth model when no population balance model"
        << "is active" << exit(FatalError);

    return nucleationRateModel::New(*this, *this)();
}


const Foam::nucleateSizeModel& Foam::populationBalances::noPBM::nucleateSize()
const
{
    FatalErrorInFunction
        << "Invalid request of growth model when no population balance model"
        << "is active" << exit(FatalError);

    return nucleateSizeModel::New(*this, *this)();
}


const Foam::surfaceScalarField& Foam::populationBalances::noPBM::phi()
const
{
    FatalErrorInFunction
        << "Invalid request of phi when no population balance model is active"
        << exit(FatalError);

    return surfaceScalarField::New("dummy", mesh_, 0.0);
}


const Foam::PtrList<Foam::volScalarField>&
Foam::populationBalances::noPBM::moments() const
{
    FatalErrorInFunction
        << "Invalid request of moments when no population balance model is "
        << "active"
        << exit(FatalError);

    return moments_;
}


Foam::PtrList<Foam::volScalarField>&
Foam::populationBalances::noPBM::moments()
{
    FatalErrorInFunction
        << "Invalid request of moments when no population balance model is "
        << "active"
        << exit(FatalError);

    return moments_;
}


const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::populationBalances::noPBM::precRate() const
{
    return precRate_;
}


const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::populationBalances::noPBM::superSat() const
{
    FatalErrorInFunction
        << "Invalid request of supersaturation when no population balance "
        << "model is active"
        << exit(FatalError);

    return precRate_;
}


const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::populationBalances::noPBM::cationConcRatios() const
{
    FatalErrorInFunction
        << "Invalid request of cationConcRatios when no population balance "
        << "model is active"
        << exit(FatalError);

    return cationConcRatios_;
}


const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::populationBalances::noPBM::sourceCorrCoeff() const
{
    FatalErrorInFunction
        << "Invalid request of sourceCorrCoeff when no population balance "
        << "model is active"
        << exit(FatalError);

    return precRate_;
}


const uint8_t& Foam::populationBalances::noPBM::numOfMoments() const
{
    return numOfMoments_;
}


const uint8_t& Foam::populationBalances::noPBM::numOfNodes() const
{
    return numOfNodes_;
}


const uint8_t& Foam::populationBalances::noPBM::numOfCoord() const
{
    return numOfCoord_;
}


const Foam::List<Foam::scalar>& Foam::populationBalances::noPBM::dSmall() const
{
    FatalErrorInFunction
        << "Invalid request of dSmall when no population balance model is "
        << "active"
        << exit(FatalError);
    
    return dSmall_;
}


Foam::scalar Foam::populationBalances::noPBM::crystalRho() const
{
    FatalErrorInFunction
        << "Invalid request of crystalRho when no population balance model is "
        << "active"
        << exit(FatalError);
    
    return 0.0;
}


Foam::scalar Foam::populationBalances::noPBM::crystalMW() const
{
    FatalErrorInFunction
        << "Invalid request of crystalMW when no population balance model is "
        << "active"
        << exit(FatalError);
    
    return 0.0;
}


// ************************************************************************* //
