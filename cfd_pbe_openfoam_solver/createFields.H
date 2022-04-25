Foam::Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
bool returnValue(setRefCell(p, piso.dict(), pRefCell, pRefValue));
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

incompressible::momentumTransportModel& turb_ref = turbulence();

solutionNMC solution_nmc(mesh, phi, turb_ref);

autoPtr<populationBalance> pb
(
    populationBalance::New(mesh, turb_ref, solution_nmc)
);

populationBalance& pb_ref = pb();

odeSolver ode(mesh, turb_ref, pb_ref, solution_nmc);

PtrList<volScalarField>& moments = pb_ref.moments();

PtrList<volScalarField>& totalMetalConcs = solution_nmc.totalMetalConcs();
// const volScalarField& totalNH3 = solution_nmc.totalNH3();
// const volScalarField& totalNa = solution_nmc.totalNa();
// const volScalarField& totalSO4 = solution_nmc.totalSO4();

// const PtrList<DimensionedField<scalar, volMesh>>& metalCationConcs =
//     solution_nmc.metalCationConcs();
// const DimensionedField<scalar, volMesh>& NH3 = solution_nmc.NH3();
// const DimensionedField<scalar, volMesh>& OH = solution_nmc.OH();

const Field<scalar>& activeCellsField = solution_nmc.activeCellsField();

#include "createMRF.H"
#include "createFvOptions.H"
