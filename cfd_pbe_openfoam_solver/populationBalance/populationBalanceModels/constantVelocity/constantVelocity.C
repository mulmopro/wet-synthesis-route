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

#include "constantVelocity.H"
#include "solutionNMC.H"
#include "quadratureMethod.H"
#include "growthModel.H"
#include "nucleationRateModel.H"
#include "nucleateSizeModel.H"
#include "aggregationList.H"
#include "breakageList.H"
#include "fvcFlux.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalances
{
    defineTypeNameAndDebug(constantVelocity, 0);
    addToRunTimeSelectionTable(populationBalance, constantVelocity, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalances::constantVelocity::constantVelocity
(
    const fvMesh& mesh,
    const incompressible::momentumTransportModel& turbulence,
    const solutionNMC& solution
)
:
    populationBalance(mesh, turbulence, solution),

    phi_(turbulence.phi()),

    superSat_(solution.superSat()),

    cationConcRatios_(solution.cationConcRatios()),

    alphaMin_("alphaMinPB", dimless, lookup("alphaMin")),

    dSmall_(lookup("dSmall")),

    numOfNodes_(readInt(lookup("numOfNodes"))),

    numOfMoments_(2*numOfNodes_),

    numOfCoord_(readInt(subDict("coordinates").lookup("numOfCoordinates"))),

    writeSummaryInterval_(lookupOrDefault("writeSummaryInterval", 1)),

    crystalRho_
    (
        "crystalRho",
        dimDensity,
        subDict("crystalProperties").lookup("crystalDensity")
    ),

    crystalMW_
    (
        "crystalMW",
        dimMass/dimMoles,
        subDict("crystalProperties").lookup("crystalMW")
    ),

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

    sourceCorrCoeff_
    (
        IOobject
        (
            "sourceCorrCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    )
{
    if (numOfNodes_ < 2)
    {
        FatalErrorInFunction
            << "Number of nodes, " << Foam::name(numOfNodes_)
            << ", is not valid."
            << endl << "Number of nodes has to be equal ro greater than 2."
            << endl << endl
            << exit(FatalError);
    }

    if (dSmall_.size() != numOfNodes_)
    {
        FatalErrorInFunction
            << "Number of small diameters, " << dSmall_ //Foam::name(dSmall_.size())
            << ", is not consistent with the number of nodes." << endl
            << "Please check the list 'dSmall' in the file 'pbPrperties'."
            << endl << endl
            << exit(FatalError);
    }
    
    for (label count=0; count < numOfMoments_; count++)
    {
        word momentName("M" + Foam::name(count));
        Info<< "Reading field " << momentName << "\n" << endl;
        moments_.append
        (
            new volScalarField
            (
                IOobject
                (
                    momentName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        ); // Declaring the moments pointer list
    }

    forAll(moments_, momenti)
    {
        mFields_.add(moments_[momenti]);
    }

    // Growth model
    growth_.set
    (
        growthModel::New
        (
            subDict("sources").subDict("growth"),
            *this
        ).ptr()
    );

    // Nucleation rate
    nucleationRate_.set
    (
        nucleationRateModel::New
        (
            subDict("sources").subDict("nucleation"),
            *this
        ).ptr()
    );

    // Nucleate size
    nucleateSize_.set
    (
        nucleateSizeModel::New
        (
            subDict("sources").subDict("nucleation").subDict("nucleateSize"),
            *this
        ).ptr()
    );
    
    // List of aggregation models
    aggregationList_.set
    (
        new aggregationList
        (
            subDict("sources").subDict("aggregation"),
            *this
        )
    );
    
    // List of breakage models
    breakageList_.set
    (
        new breakageList
        (
            subDict("sources").subDict("breakage"),
            *this
        )
    );
    
    // Quadrature method
    quadrature_.set
    (
        quadratureMethod::New(*this).ptr()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalances::constantVelocity::~constantVelocity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::populationBalances::constantVelocity::transport_moments()
{

    // const DimensionedField<scalar, volMesh>& M3 = moments_[3].internalField();

    sourceCorrCoeff_ = pos(moments_[3].internalField()*kv_ - alphaMin_);
    forAll(moments_, momenti)
    {
        if (momenti != 3)
        {
            sourceCorrCoeff_ *= pos(moments_[momenti].internalField());
        }
    }

    Info<< "Calculating nodes and weights in "
        << gSum(sourceCorrCoeff_)
        << " number of cells"<< endl;

    // Finding nodes and weights for tagged cells
    quadrature_->nodesAndWeights(sourceCorrCoeff_.field());
    
    const PtrList<volScalarField>& nodes = quadrature_->nodes();
    const PtrList<volScalarField>& weights = quadrature_->weights();

    // Moment transport equations without growth and nucleation
    #include "momentTransportEq.H"
}


void Foam::populationBalances::constantVelocity::correct()
{
    bool writeSummary = !(mesh_.time().timeIndex() % writeSummaryInterval_);

    forAll(moments_, momenti)
    {
        volScalarField& M = moments_[momenti];

        M.correctBoundaryConditions();

        if (writeSummary)
        {
            Info<< M.name() << " = "
                << M.weightedAverage(mesh_.V()).value()
                << "  Min(" << M.name() << ") = " << min(M).value()
                << "  Max(" << M.name() << ") = " << max(M).value()
                << endl;
        }

        // M.max(0);
        // M.correctBoundaryConditions();
    }
}


bool Foam::populationBalances::constantVelocity::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        writeSummaryInterval_ = this->lookupOrDefault(
            "writeSummaryInterval", 1);

        // this->lookup("alphaMin") >> alphaMin_;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
