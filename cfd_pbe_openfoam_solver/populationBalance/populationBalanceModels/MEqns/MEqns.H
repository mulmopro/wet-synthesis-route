tmp<fv::convectionScheme<scalar>> mvConv
(
    fv::convectionScheme<scalar>::New
    (
        mesh_, mFields_, phi_, mesh_.divScheme("div(phi,M)")
    )
);

PtrList<DimensionedField<scalar, volMesh>> growthByWeights;

forAll(nodes, nodei)
{
    growthByWeights.append
    (
        new DimensionedField<scalar, volMesh>
        (
            "growthByWeight" + Foam::name(nodei),
            growth_->rate(nodes[nodei])*weights[nodei]
        )
    );
}


tmp<DimensionedField<scalar, volMesh>> tNucRate
(
    new DimensionedField<scalar, volMesh>
    (
        "nucRate",
        nucleationRate_->rate()
    )
);

DimensionedField<scalar, volMesh>& nucRate = tNucRate.ref();


tmp<DimensionedField<scalar, volMesh>> tNucSize
(
    new DimensionedField<scalar, volMesh>
    (
        "nucSize",
        nucleateSize_->size()
    )
);

DimensionedField<scalar, volMesh>& nucSize = tNucSize.ref();

// Update precRate
tmp<DimensionedField<scalar, volMesh>> tx
(
    new DimensionedField<scalar, volMesh>
    (
        IOobject
        (
            "sourceM3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    )
);

DimensionedField<scalar, volMesh>& x = tx.ref();

forAll(nodes, nodei)
{
    x += growthByWeights[nodei]
       * Foam::pow(nodes[nodei].internalField(), 2);
}

x *= 3.0;

x += nucRate * Foam::pow(nucSize, 3);

precRate_ = kv_ * x * crystalRho_ / crystalMW_;

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
        forAll(nodes, nodei)
        {
            sourceM += growthByWeights[nodei]
                * Foam::pow(nodes[nodei].internalField(), momenti - 1);
        }

        sourceM *= momenti;

        sourceM +=
            nucRate * Foam::pow(nucSize, momenti)
          + aggregationList_->source(nodes, weights, momenti)
          + breakageList_->source(nodes, weights, momenti);
    }
    else
    {
        sourceM += x;
    }

    // Definition of the moment equation corresponding to the current moment
    tmp<fvScalarMatrix> tMiEqn
    (
        fvm::ddt(M)
    //   + fvm::div(phi_, M)
    //   +  mvConv->fvcDiv(phi_, M)
      + mvConv->fvmDiv(phi_, M)
      ==
        sourceM
    );

    fvScalarMatrix& MiEqn = tMiEqn.ref();

    // Relaxing the moment equation
    MiEqn.relax();
    
    // Solve the equation of the current moment
    MiEqn.solve(M.name());

    M.correctBoundaryConditions();

    Info<< M.name() << " = "
        << M.weightedAverage(mesh_.V()).value()
        << "  Min(" << M.name() << ") = " << min(M).value()
        << "  Max(" << M.name() << ") = " << max(M).value()
        << endl;
    
    // M.max(0);
    // M.correctBoundaryConditions();
}

tx.clear();
