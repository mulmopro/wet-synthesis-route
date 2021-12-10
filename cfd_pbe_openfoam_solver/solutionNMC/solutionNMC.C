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

#include "solutionNMC.H"
#include "linearSolver.H"
#include "CEqn.H"

#include "UserData.H"

#include "cellSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solutionNMC::equilibriaEqs
(
    const List<scalar>& pConcs, const List<scalar>& totalConcs,
    scalar cationTotalConc, double* negf, double* J
)
{
    label i, j;
    label countComplex = 0;

    scalar conc_OH = Foam::pow(10, -1.0*pConcs[indexOH_]);
    scalar conc_NH3 = Foam::pow(10, -1.0*pConcs[indexNH3_]);
    scalar conc_NH4 = Kb_NH3_*conc_NH3 / conc_OH;

    scalar f_NH3 = 0.0;
    scalar J_NH3 = 0.0;
    for (i=0; i<nMetals_; i++)
    {
        label numOfComplex = numOfComplexes_[i];

        scalar conc_i = Foam::pow(10, -1.0*pConcs[i]);

        scalar f_i = 0.0;
        scalar f_NH3_i = 0.0;
        for (j=1; j<=numOfComplex; j++)
        {
            scalar temp = kn_NMC_[countComplex]*conc_i*Foam::pow(conc_NH3, j);
            f_i += temp;
            f_NH3_i += j*temp;
            J_NH3 += Foam::pow(j, 2)*temp;

            countComplex ++;
        }

        negf[i] = conc_i + f_i - totalConcs[i];

        J[i*nComps_ + i] = ln10_*(conc_i + f_i);

        /* The index is changed with respect to the python code
        to respect the order required by the dgesv routine of LAPACK */
        J[indexNH3_*nComps_ + i] = ln10_*f_NH3_i; 
        J[i*nComps_ + indexNH3_] = ln10_*f_i;

        f_NH3 += f_NH3_i;
    }

    negf[indexNH3_] = conc_NH3 + conc_NH4 + f_NH3 - totalConcs[indexNH3_];

    J[indexNH3_*nComps_ + indexNH3_] = ln10_*(conc_NH3 + conc_NH4 + J_NH3);
    J[indexOH_*nComps_ + indexNH3_] = -1.0*ln10_*conc_NH4;

    negf[indexOH_] = conc_OH
        - totalConcs[indexOH_] - 2*cationTotalConc - conc_NH4 - kw_/conc_OH;

    J[indexOH_*nComps_ + indexOH_] = ln10_*(conc_NH4 + kw_/conc_OH + conc_OH);
    J[indexNH3_*nComps_ + indexOH_] = -1.0*ln10_*conc_NH4;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionNMC::solutionNMC
(
    const fvMesh& mesh,
    const surfaceScalarField& phi,
    const incompressible::momentumTransportModel& turbulence
)
:
    IOdictionary
    (
        IOobject
        (
            "solutionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phi_(phi),

    turbulence_(turbulence),

    metalNames_(lookup("metals")),

    nMetals_(metalNames_.size()),

    nComps_(nMetals_ + 2),

    indexNH3_(nMetals_),

    indexOH_(nComps_ - 1),

    maxIter_(readInt(subDict("newtonRaphson").lookup("maxIter"))),

    tol_(readScalar(subDict("newtonRaphson").lookup("tol"))),

    smallConc_(readScalar(lookup("smallConc"))),

    effectiveConc_(readScalar(lookup("effectiveConc"))),

    writeSummaryInterval_(lookupOrDefault("writeSummaryInterval", 1)),

    turbSc_
    (
        dimensionedScalar
        (
            "turbSc_", dimless, lookup("turbulentSchmidt")
        )
    ),

    activeCells_
    (
        IOobject
        (
            "activeCells",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("one", dimless, 1.0),
        false
    ),

    activeCellsField_(activeCells_.field()),

    NH3_
    (
        IOobject
        (
            "NH3",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMoles/dimVol, 0.0)
    ),

    OH_
    (
        IOobject
        (
            "OH",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMoles/dimVol, 0.0)
    ),

    superSat_
    (
        IOobject
        (
            "supersaturation",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    pH_
    (
        IOobject
        (
            "pH",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    
    totalNH3_
    (
        IOobject
        (
            "total_NH3",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    inertCharges_
    (
        IOobject
        (
            "inertCharges",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    for (label count=0; count < nMetals_; count++)
    {
        word metalName(metalNames_[count]);

        Info<< "Reading field total_" << metalName << "\n" << endl;
        totalMetalConcs_.append
        (
            new volScalarField
            (
                IOobject
                (
                    "total_" + metalName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        metalCationConcs_.append
        (
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    metalName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimConc, 0.0)
            )
        );

        cationConcRatios_.append
        (
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "ratio_" + metalName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    List<scalar> DList(lookup("molecularDiffusivity"));

    if (DList.size() == nComps_)
    {
        for (label count=0; count < nComps_; count++)
        {
            D_.append
            (
                new dimensionedScalar
                (
                    "D_" + Foam::name(count),
                    dimArea/dimTime,
                    DList[count]
                )
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "The number of molecular diffusivities does not match "
            << "the number of transported scalars (" << nComps_ << ")"
            << exit(FatalError);
    }

    List<word> inactiveRegionNames
    (
        lookupOrDefault("chemicallyInactiveRegions", wordList())
    );

    if (!inactiveRegionNames.empty())
    {
        cellSet inactiveCellSet(mesh_, inactiveRegionNames[0]);

        forAll(inactiveRegionNames, i)
        {
            if (i > 0)
            {
                cellSet region(mesh_, inactiveRegionNames[i]);
                inactiveCellSet.addSet(region);
            }
        }

        Info<< inactiveCellSet.size()<< " cells are excluded from equilibrium "
            << "calculations"<< endl<< endl;

        forAll(activeCells_, celli)
        {
            if (inactiveCellSet.found(celli))
            {
                activeCells_[celli] = 0.0;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionNMC::~solutionNMC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solutionNMC::update(const realtype *y, UserData *aux_data,
    realtype *cationConcRatios)
{
    double negf[nComps_];
    double J[nComps_*nComps_];

    label i, j;
    int n = nComps_, nrhs = 1, lda = nComps_, ldb = nComps_, info;
    int ipiv[nComps_];

    scalar error, equilConc, totalConc, cationTotalConc, cationConcRatio,
        cationConc, conc_OH, conc_Na, k_sp_NMC, powConcs_NMC;

    Switch validConc(true);

    List<scalar> pConcs(nComps_);
    List<scalar> totalConcs(nComps_);

    for (i=0; i<nMetals_; i++)
    {
        totalConc = y[i];

        if (totalConc < effectiveConc_)
        {
            validConc = false;
            break;
        }

        totalConcs[i] = totalConc;
    }

    if (validConc)
    {
        totalConcs[indexNH3_] = aux_data->totalNH3;
        totalConcs[indexOH_] = aux_data->inertCharges;

        if (totalConcs[indexNH3_] > effectiveConc_)
        {
            cationTotalConc = 0.0;
            for (i=0; i<nMetals_; i++)
            {
                equilConc = aux_data->metalCationConcs[i];
                totalConc = totalConcs[i];
                if (equilConc > 0.0 && equilConc < totalConc)
                {
                    pConcs[i] = -1.0*Foam::log10(equilConc);
                }
                else
                {
                    pConcs[i] = -1.0*Foam::log10(totalConc);
                }
                cationTotalConc += totalConc;
            }

            equilConc = aux_data->concOH;
            conc_Na = totalConcs[indexOH_] + 2*cationTotalConc;
            if (equilConc > 1e-7)
            {
                pConcs[indexOH_] = -1.0*Foam::log10(equilConc);
            }
            else if (conc_Na > 0.001)
            {
                pConcs[indexOH_] = -1.0*Foam::log10(conc_Na);
            }
            else
            {
                pConcs[indexOH_] = -1.0*Foam::log10(0.001);
            }

            equilConc = aux_data->concNH3;
            totalConc = totalConcs[indexNH3_];
            if (equilConc > effectiveConc_ && equilConc < totalConc)
            {
                pConcs[indexNH3_] = -1.0*Foam::log10(equilConc);
            }
            else
            {
                pConcs[indexNH3_] = -1.0*Foam::log10(totalConc);
            }

            for (i=0; i<nComps_; i++)
            {
                negf[i] = 0.0;
            }

            for (i=0; i<nComps_*nComps_; i++)
            {
                J[i] = 0.0;
            }
            
            for(i=0; i<maxIter_; i++)
            {
                equilibriaEqs(pConcs, totalConcs, cationTotalConc,
                    negf, J);

                dgesv_(&n, &nrhs, J, &lda, ipiv, negf, &ldb, &info);

                if (info == 0)
                {
                    error = 0.0;
                    for (j=0; j<nComps_; j++)
                    {
                        error += mag(negf[j]);
                    }

                    if (error < tol_)
                    {
                        break;
                    }
                    else
                    {
                        for (j=0; j<nComps_; j++)
                        {
                            pConcs[j] += negf[j];
                        }
                    }

                }
                else if (info > 0) /* Check for the exact singularity */
                {
                    Info<< "The diagonal element of the triangular"
                        << "factor of A," << endl
                        << "U(" << info << "," << info << ") is zero, "
                        << "so that A is singular;" << endl
                        << "the solution could not be computed."
                        << endl << endl
                        << exit(FatalError);
                }
                else
                {
                    Info<< "Equilibrium calculation failed for cell "
                        << aux_data->cell_id << "with the total concentrations:"
                        << endl << totalConcs << endl << endl
                        << exit(FatalError);
                }
            }

            k_sp_NMC = 1.0;
            powConcs_NMC = 1.0;
            for (j=0; j<nMetals_; j++)
            {
                cationConc = Foam::pow(10, -1.0*pConcs[j]);
                aux_data->metalCationConcs[j] = cationConc;

                cationConcRatio = totalConcs[j] / cationTotalConc;
                cationConcRatios[j] = cationConcRatio;

                k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);
                powConcs_NMC *= Foam::pow(cationConc, cationConcRatio);
            }

            aux_data->concNH3 = Foam::pow(10, -1.0*pConcs[indexNH3_]);

            conc_OH = Foam::pow(10, -1.0*pConcs[indexOH_]);
            aux_data->concOH = conc_OH;

            aux_data->superSat =
                Foam::pow
                (
                    powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                    1.0/3.0
                );
        }
        else // totalConcs[indexNH3_] <= effectiveConc_
        {
            aux_data->concNH3 = totalConcs[indexNH3_];

            cationTotalConc = 0.0;
            for (i=0; i<nMetals_; i++)
            {
                cationTotalConc += totalConcs[i];
            }

            conc_OH = totalConcs[indexOH_] + 2*cationTotalConc;

            if (conc_OH > 1e-7)
            {
                k_sp_NMC = 1.0;
                powConcs_NMC = 1.0;
                for (j=0; j<nMetals_; j++)
                {
                    cationConcRatio = totalConcs[j] / cationTotalConc;
                    cationConcRatios[j] = cationConcRatio;

                    k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);
                    powConcs_NMC *=
                        Foam::pow(totalConcs[j], cationConcRatio);
                }

                aux_data->superSat =
                    Foam::pow
                    (
                        powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                        1.0/3.0
                    );

                aux_data->concOH = conc_OH;
            }
            else
            {
                for (j=0; j<nMetals_; j++)
                {
                    cationConcRatios[j] = totalConcs[j] / cationTotalConc;
                }

                aux_data->superSat = 0.0;
                aux_data->concOH = 1e-7;
            }
        }
    }
    else // validConc = False
    {
        for (j=0; j<nMetals_; j++)
        {
            aux_data->metalCationConcs[j] = 0.0;
            cationConcRatios[j] = 0.0;
        }
        aux_data->concNH3 = 0.0;
        aux_data->concOH = 0.0;
        aux_data->superSat = 0.0;
    }
}

void Foam::solutionNMC::update()
{
    double negf[nComps_];
    double J[nComps_*nComps_];

    label i, j;
    int n = nComps_, nrhs = 1, lda = nComps_, ldb = nComps_, info;
    int ipiv[nComps_];

    scalar error, equilConc, totalConc, cationTotalConc, cationConcRatio,
        cationConc, conc_OH, conc_Na, k_sp_NMC, powConcs_NMC;

    Switch validConc;

    List<scalar> pConcs(nComps_);
    List<scalar> totalConcs(nComps_);

    forAll(superSat_, celli)
    {
        if (activeCellsField_[celli] > 0)
        {
            validConc = true;
            cationTotalConc = 0.0;
            for (i=0; i<nMetals_; i++)
            {
                totalConc = totalMetalConcs_[i][celli];

                if (totalConc < effectiveConc_)
                {
                    validConc = false;
                    // break;
                }

                totalConcs[i] = totalConc;
                cationTotalConc += totalConc;
            }

            totalConcs[indexNH3_] = totalNH3_[celli];
            totalConcs[indexOH_] = inertCharges_[celli];

            if (validConc)
            {
                // totalConcs[indexNH3_] = totalNH3_[celli];
                // totalConcs[indexOH_] = inertCharges_[celli];

                if (totalConcs[indexNH3_] > effectiveConc_)
                {
                    // cationTotalConc = 0.0;
                    for (i=0; i<nMetals_; i++)
                    {
                        equilConc = metalCationConcs_[i][celli];
                        totalConc = totalConcs[i];
                        if (equilConc > 0.0 && equilConc < totalConc)
                        {
                            pConcs[i] = -1.0*Foam::log10(equilConc);
                        }
                        else
                        {
                            pConcs[i] = -1.0*Foam::log10(totalConc);
                        }
                        // cationTotalConc += totalConc;
                    }

                    equilConc = OH_[celli];
                    conc_Na = totalConcs[indexOH_] + 2*cationTotalConc;
                    if (equilConc > 1e-7)
                    {
                        pConcs[indexOH_] = -1.0*Foam::log10(equilConc);
                    }
                    else if (conc_Na > 0.001)
                    {
                        pConcs[indexOH_] = -1.0*Foam::log10(conc_Na);
                    }
                    else
                    {
                        pConcs[indexOH_] = -1.0*Foam::log10(0.001);
                    }

                    equilConc = NH3_[celli];
                    totalConc = totalConcs[indexNH3_];
                    if (equilConc > effectiveConc_ && equilConc < totalConc)
                    {
                        pConcs[indexNH3_] = -1.0*Foam::log10(equilConc);
                    }
                    else
                    {
                        pConcs[indexNH3_] = -1.0*Foam::log10(totalConc);
                    }

                    for (i=0; i<nComps_; i++)
                    {
                        negf[i] = 0.0;
                    }

                    for (i=0; i<nComps_*nComps_; i++)
                    {
                        J[i] = 0.0;
                    }
                    
                    for(i=0; i<maxIter_; i++)
                    {
                        equilibriaEqs(pConcs, totalConcs, cationTotalConc,
                            negf, J);

                        dgesv_(&n, &nrhs, J, &lda, ipiv, negf, &ldb, &info);

                        if (info == 0)
                        {
                            error = 0.0;
                            for (j=0; j<nComps_; j++)
                            {
                                error += mag(negf[j]);
                            }

                            if (error < tol_)
                            {
                                break;
                            }
                            else
                            {
                                for (j=0; j<nComps_; j++)
                                {
                                    pConcs[j] += negf[j];
                                }
                            }

                        }
                        else if (info > 0) /* Check for the exact singularity */
                        {
                            Info<< "The diagonal element of the triangular"
                                << "factor of A," << endl
                                << "U(" << info << "," << info << ") is zero, "
                                << "so that A is singular;" << endl
                                << "the solution could not be computed."
                                << endl << endl
                                << exit(FatalError);
                        }
                        else
                        {
                            Info<< "Equilibrium calculation failed for cell "
                                << celli << "with the total concentrations:"
                                << endl << totalConcs << endl << endl
                                << exit(FatalError);
                        }
                    }

                    k_sp_NMC = 1.0;
                    powConcs_NMC = 1.0;
                    for (j=0; j<nMetals_; j++)
                    {
                        cationConc = Foam::pow(10, -1.0*pConcs[j]);
                        metalCationConcs_[j][celli] = cationConc;

                        cationConcRatio = totalConcs[j] / cationTotalConc;
                        cationConcRatios_[j][celli] = cationConcRatio;
                        
                        k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);
                        powConcs_NMC *= Foam::pow(cationConc, cationConcRatio);
                    }

                    NH3_[celli] = Foam::pow(10, -1.0*pConcs[indexNH3_]);

                    conc_OH = Foam::pow(10, -1.0*pConcs[indexOH_]);
                    OH_[celli] = conc_OH;

                    pH_[celli] = 14 - pConcs[indexOH_];

                    superSat_[celli] =
                        Foam::pow
                        (
                            powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                            1.0/3.0
                        );
                }
                else // totalConcs[indexNH3_] <= effectiveConc_
                {
                    NH3_[celli] = totalConcs[indexNH3_];

                    // cationTotalConc = 0.0;
                    // for (i=0; i<nMetals_; i++)
                    // {
                    //     cationTotalConc += totalConcs[i];
                    // }

                    conc_OH = totalConcs[indexOH_] + 2*cationTotalConc;

                    if (conc_OH > 1e-7)
                    {
                        k_sp_NMC = 1.0;
                        powConcs_NMC = 1.0;
                        for (j=0; j<nMetals_; j++)
                        {
                            cationConcRatio = totalConcs[j] / cationTotalConc;
                            cationConcRatios_[j][celli] = cationConcRatio;

                            k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);
                            powConcs_NMC *=
                                Foam::pow(totalConcs[j], cationConcRatio);
                        }

                        superSat_[celli] =
                            Foam::pow
                            (
                                powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                                1.0/3.0
                            );

                        OH_[celli] = conc_OH;
                        pH_[celli] = 14 + Foam::log10(conc_OH);
                    }
                    else
                    {
                        for (j=0; j<nMetals_; j++)
                        {
                            cationConcRatios_[j][celli] = 
                                totalConcs[j] / cationTotalConc;
                        }

                        superSat_[celli] = 0.0;
                        OH_[celli] = 1e-7;
                        pH_[celli] = 7.0;
                    }
                }
            }
            else // validConc = False
            {
                for (j=0; j<nMetals_; j++)
                {
                    metalCationConcs_[j][celli] = 0.0;
                    cationConcRatios_[j][celli] = 0.0;
                }
                NH3_[celli] = 0.0;

                superSat_[celli] = 1.0;

                conc_OH =
                    (
                        Foam::sqrt
                        (
                            4
                            *(
                                kw_
                                +   Foam::max
                                    (
                                        Kb_NH3_*totalConcs[indexNH3_], 0.0
                                    )
                            )
                            + Foam::sqr(totalConcs[indexOH_])
                        )
                        + totalConcs[indexOH_]
                    ) / 2.0;

                OH_[celli] = 0.0;
                pH_[celli] = 14 + Foam::log10(conc_OH);
            }
        }
        else // cell is not in the chemically active zone
        {
            for (j=0; j<nMetals_; j++)
            {
                metalCationConcs_[j][celli] = 0.0;
                cationConcRatios_[j][celli] = 0.0;
            }
            NH3_[celli] = 0.0;
            OH_[celli] = 0.0;
            superSat_[celli] = 0.0;
            pH_[celli] = 0.0;
        }
    }
}


void Foam::solutionNMC::transport_species()
{
    volScalarField Dturb = turbulence_.nut() / turbSc_;

    // A temporary object to define the concentration transport equations
    tmp<fvScalarMatrix> tCiEqn;

    // Loop over metals
    forAll(totalMetalConcs_, metali)
    {
        // Reference to the metal concentration for which transport eq. will be
        // defined
        volScalarField& metal = totalMetalConcs_[metali];

        volScalarField DEff = Dturb + D_[metali];

        makeConcTransportEq(metal);
    }

    volScalarField DEff = Dturb + D_[indexNH3_];

    makeConcTransportEq(totalNH3_);

    DEff = Dturb + D_[indexOH_];

    makeConcTransportEq(inertCharges_);

}


void Foam::solutionNMC::correct()
{
    bool writeSummary = !(mesh_.time().timeIndex() % writeSummaryInterval_);

    // Loop over metals
    forAll(totalMetalConcs_, metali)
    {
        // Reference to the metal concentration
        volScalarField& metal = totalMetalConcs_[metali];

        correctSpecies(metal, writeSummary);
    }

    correctSpecies(totalNH3_, writeSummary);
    correctSpecies(inertCharges_, writeSummary);
}


bool Foam::solutionNMC::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        writeSummaryInterval_ = this->lookupOrDefault(
            "writeSummaryInterval", 1);

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
