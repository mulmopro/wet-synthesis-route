#ifndef makeConcTransportEq_H
#define makeConcTransportEq_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeConcTransportEq(Specie)                                           \
                                                                              \
    tCiEqn =                                                                  \
    (                                                                         \
        fvm::ddt(Specie)                                                      \
      + fvm::div(phi_, Specie)                                                \
      - fvm::laplacian(DEff, Specie)                                          \
    );                                                                        \
                                                                              \
    tCiEqn.ref().relax();                                                     \
                                                                              \
    tCiEqn.ref().solve(Specie.name());                                        \
                                                                              \
    /* Specie.correctBoundaryConditions();                                       \
                                                                              \
    Info<< Specie.name() << " = "                                             \
        << Specie.weightedAverage(mesh_.V()).value()                          \
        << "  Min(" << Specie.name() << ") = " << min(Specie).value()         \
        << "  Max(" << Specie.name() << ") = " << max(Specie).value()         \
        << endl;                                                              \ */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


#ifndef correctSpecies_H
#define correctSpecies_H

// #include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define correctSpecies(Specie)                                                \
                                                                              \
    Specie.correctBoundaryConditions();                                       \
                                                                              \
    Info<< Specie.name() << " = "                                             \
        << Specie.weightedAverage(mesh_.V()).value()                          \
        << "  Min(" << Specie.name() << ") = " << min(Specie).value()         \
        << "  Max(" << Specie.name() << ") = " << max(Specie).value()         \
        << endl;                                                              \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
