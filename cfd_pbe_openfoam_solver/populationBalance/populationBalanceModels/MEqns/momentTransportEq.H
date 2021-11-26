tmp<fv::convectionScheme<scalar>> mvConv
(
    fv::convectionScheme<scalar>::New
    (
        mesh_, mFields_, phi_, mesh_.divScheme("div(phi,M)")
    )
);

volScalarField Dturb = turbulence_.nut();

// Loop over moments to define their corresponding moment equations
forAll(moments_, momenti)
{
    // Reference to the moment for which transport eq. is going to be defined
    volScalarField& M = moments_[momenti];

    // Calculate the source of the current moment 
    DimensionedField<scalar, volMesh> sourceM
    (
        IOobject
        (
            "sourceM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", M.dimensions()/dimTime, 0.0)
    );

    if (momenti != 3)
    {
        sourceM +=
            aggregationList_->source(nodes, weights, momenti)
          + breakageList_->source(nodes, weights, momenti);
    }

    // Definition of the moment equation corresponding to the current moment
    tmp<fvScalarMatrix> tMiEqn
    (
        fvm::ddt(M)
    //   + fvm::div(phi_, M)
    //   +  mvConv->fvcDiv(phi_, M)
      + mvConv->fvmDiv(phi_, M)
      - fvm::laplacian(Dturb, M)
      ==
        sourceM
    );

    fvScalarMatrix& MiEqn = tMiEqn.ref();

    // Relaxing the moment equation
    MiEqn.relax();
    
    // Solve the equation of the current moment
    MiEqn.solve(M.name());

    // M.correctBoundaryConditions();

    // Info<< M.name() << " = "
    //     << M.weightedAverage(mesh_.V()).value()
    //     << "  Min(" << M.name() << ") = " << min(M).value()
    //     << "  Max(" << M.name() << ") = " << max(M).value()
    //     << endl;

    M.max(0);
    // M.correctBoundaryConditions();
}
