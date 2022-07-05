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
#include "activityCoeffModel.H"
#include "cellSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solutionNMC::equilibriaEqs
(
    const List<scalar>& pConcs, const List<scalar>& totalConcs,
    scalar cationTotalConc, double* negf, double* J
) const
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

    negf[indexOH_] = conc_OH + 2*totalConcs[indexSO4_]
        - totalConcs[indexNa_] - 2*cationTotalConc - conc_NH4 - kw_/conc_OH;

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

    nTotalConc_(nMetals_ + 3),

    nComps_(nMetals_ + 2),

    indexNH3_(nMetals_),

    indexOH_(nComps_ - 1),

    indexNa_(indexNH3_ + 1),

    indexSO4_(indexNa_ + 1),

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

    totalNa_
    (
        IOobject
        (
            "total_Na",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    totalSO4_
    (
        IOobject
        (
            "total_SO4",
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

    if (DList.size() == nTotalConc_)
    {
        for (label count=0; count < nTotalConc_; count++)
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
            << "the number of transported scalars (" << nTotalConc_ << ")"
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

    activityModel_.set
    (
        activityCoeffModel::New
        (
            subDict("activity"),
            *this
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionNMC::~solutionNMC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solutionNMC::solve(const realtype *y, UserData *aux_data,
    realtype *cationConcRatios) const
{
    label nJ = nComps_*nComps_;

    double negf[nComps_];
    double J[nJ];

    label i, j, k;
    int n = nComps_, nrhs = 1, lda = nComps_, ldb = nComps_, info;
    int ipiv[nComps_];

    scalar error, equilConc, totalConc, cationConcRatio, conc_NH3, conc_OH,
        conc_NH4, conc_H, k_sp_NMC, powConcs_NMC;

    List<scalar> totalConcs(nTotalConc_);
    List<scalar> pConcs(nComps_);

    scalar cationTotalConc(0.0);
    for (i=0; i<nMetals_; i++)
    {
        totalConc = y[i];

        totalConcs[i] = totalConc;
        cationTotalConc += totalConc;
    }

    scalar totalConc_NH3(aux_data->totalNH3);
    scalar totalConc_Na(aux_data->totalNa);
    scalar totalConc_SO4(aux_data->totalSO4);

    if (totalConc_NH3 > effectiveConc_)
    {
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
        }

        equilConc = aux_data->concNH3;
        if (equilConc > effectiveConc_ && equilConc < totalConc_NH3)
        {
            pConcs[indexNH3_] = -1.0*Foam::log10(equilConc);
        }
        else
        {
            pConcs[indexNH3_] = -1.0*Foam::log10(totalConc_NH3);
        }

        equilConc = aux_data->concOH;
        if (equilConc > 1e-7)
        {
            pConcs[indexOH_] = -1.0*Foam::log10(equilConc);
        }
        else if (totalConc_Na > 0.001)
        {
            pConcs[indexOH_] = -1.0*Foam::log10(totalConc_Na);
        }
        else
        {
            pConcs[indexOH_] = -1.0*Foam::log10(0.001);
        }

        for (i=0; i<nComps_; i++)
        {
            negf[i] = 0.0;
        }

        for (i=0; i<nJ; i++)
        {
            J[i] = 0.0;
        }

        totalConcs[indexNH3_] = totalConc_NH3;
        totalConcs[indexNa_] = totalConc_Na;
        totalConcs[indexSO4_] = totalConc_SO4;

        bool reset_guess = false;
        for (i=0; i<maxIter_; i++)
        {
            for (k=0; k<nComps_; k++)
            {
                if (pConcs[k] > 20 || pConcs[k] < -2)
                {
                    reset_guess = true;
                }
            }

            if (reset_guess)
            {
                reset_guess = false;

                if (totalConc_NH3 > 10*cationTotalConc)
                {
                    for (k=0; k<nMetals_; k++)
                    {
                        pConcs[k] = -1.0*Foam::log10(totalConcs[k]/100);
                    }

                    pConcs[indexNH3_] = -1.0*Foam::log10(
                        totalConc_NH3 - cationTotalConc);
                }
                else if (totalConc_NH3 > cationTotalConc)
                {
                    for (k=0; k<nMetals_; k++)
                    {
                        pConcs[k] = -1.0*Foam::log10(totalConcs[k]/10);
                    }

                    pConcs[indexNH3_] = -1.0*Foam::log10(
                        totalConc_NH3 - cationTotalConc / 2.0);
                }
                else
                {
                    for (k=0; k<nMetals_; k++)
                    {
                        pConcs[k] = -1.0*Foam::log10(totalConcs[k]);
                    }

                    pConcs[indexNH3_] = -1.0*Foam::log10(totalConc_NH3);
                }

                conc_OH = totalConc_Na + 2*(cationTotalConc - totalConc_SO4);
                if (conc_OH > 1e-7)
                {
                    pConcs[indexOH_] = -1.0*Foam::log10(conc_OH);
                }
                else
                {
                    pConcs[indexOH_] = -1.0*Foam::log10(0.0001);
                }
            }

            equilibriaEqs(pConcs, totalConcs, cationTotalConc, negf, J);

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

        conc_NH3 = Foam::pow(10, -1.0*pConcs[indexNH3_]);
        conc_OH = Foam::pow(10, -1.0*pConcs[indexOH_]);

        conc_NH4 = Kb_NH3_*conc_NH3 / conc_OH;
        conc_H = kw_ / conc_OH;

        List<scalar> cationMolalConc(nTotalConc_);
        for (j=0; j<nMetals_; j++)
        {
            cationMolalConc[j] = totalConcs[j];
        }
        cationMolalConc[nMetals_] = totalConc_Na;
        cationMolalConc[nMetals_ + 1] = conc_NH4;
        cationMolalConc[nMetals_ + 2] = conc_H;

        List<scalar> anionMolalConc({totalConc_SO4, conc_OH});

        scalar I_s = activityModel_->ionic_strength(
            cationMolalConc, anionMolalConc);

        k_sp_NMC = 1.0;
        powConcs_NMC = 1.0;
        for (j=0; j<nMetals_; j++)
        {
            equilConc = Foam::pow(10, -1.0*pConcs[j]);
            // aux_data->metalCationConcs[j] = equilConc;

            cationConcRatio = totalConcs[j] / cationTotalConc;
            cationConcRatios[j] = cationConcRatio;

            k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);

            scalar gamma_ca = activityModel_->pair_activity_coeff(
                j, cationMolalConc, anionMolalConc, I_s);

            powConcs_NMC *=
                Foam::pow
                (
                    equilConc * Foam::pow(gamma_ca, 3),
                    cationConcRatio
                );
        }

        // aux_data->concNH3 = conc_NH3;
        // aux_data->concOH = conc_OH;

        aux_data->superSat =
            Foam::pow
            (
                powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                1.0/3.0
            );
    }
    else // totalConc_NH3 <= effectiveConc_
    {
        // aux_data->concNH3 = totalConc_NH3;

        conc_OH = totalConc_Na + 2*(cationTotalConc - totalConc_SO4);

        if (conc_OH > 1e-7)
        {
            conc_NH4 = 0.0;
            conc_H = kw_ / conc_OH;

            List<scalar> cationMolalConc(nTotalConc_);
            for (j=0; j<nMetals_; j++)
            {
                cationMolalConc[j] = totalConcs[j];
            }
            cationMolalConc[nMetals_] = totalConc_Na;
            cationMolalConc[nMetals_ + 1] = conc_NH4;
            cationMolalConc[nMetals_ + 2] = conc_H;

            List<scalar> anionMolalConc({totalConc_SO4, conc_OH});

            scalar I_s = activityModel_->ionic_strength(
                cationMolalConc, anionMolalConc);

            k_sp_NMC = 1.0;
            powConcs_NMC = 1.0;
            for (j=0; j<nMetals_; j++)
            {
                cationConcRatio = totalConcs[j] / cationTotalConc;
                cationConcRatios[j] = cationConcRatio;

                k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);

                scalar gamma_ca = activityModel_->pair_activity_coeff(
                    j, cationMolalConc, anionMolalConc, I_s);

                powConcs_NMC *=
                    Foam::pow
                    (
                        totalConcs[j] * Foam::pow(gamma_ca, 3),
                        cationConcRatio
                    );
            }

            aux_data->superSat =
                Foam::pow
                (
                    powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                    1.0/3.0
                );

            // aux_data->concOH = conc_OH;
        }
        else
        {
            for (j=0; j<nMetals_; j++)
            {
                cationConcRatios[j] = totalConcs[j] / cationTotalConc;
            }

            aux_data->superSat = 0.0;
            // aux_data->concOH = 1e-7;
        }
    }
}

bool Foam::solutionNMC::update(label celli)
{
    label nJ = nComps_*nComps_;

    double negf[nComps_];
    double J[nJ];

    label i, j;
    int n = nComps_, nrhs = 1, lda = nComps_, ldb = nComps_, info;
    int ipiv[nComps_];

    scalar error, equilConc, totalConc, totalConc_NH3, totalConc_Na,
        totalConc_SO4, cationTotalConc, cationConcRatio, conc_NH3, conc_OH,
        conc_NH4, conc_H, charge_inerts, k_sp_NMC, powConcs_NMC;

    bool validTotalConc, validEquilConc;

    List<scalar> totalConcs(nTotalConc_);
    List<scalar> equilConcs(nComps_);
    List<scalar> pConcs(nComps_);

    validTotalConc = true;
    validEquilConc = true;

    cationTotalConc = 0.0;
    for (i=0; i<nMetals_; i++)
    {
        totalConc = totalMetalConcs_[i][celli];

        totalConcs[i] = totalConc;
        cationTotalConc += totalConc;

        if (totalConc < 0.0)
        {
            validTotalConc = false;
        }
        else
        {
            equilConc = metalCationConcs_[i][celli];

            if (equilConc > 0.0 && equilConc < totalConc)
            {
                equilConcs[i] = equilConc;
            }
            else
            {
                validEquilConc = false;
            }
        }
    }

    totalConc_NH3 = totalNH3_[celli];
    totalConc_Na = totalNa_[celli];
    totalConc_SO4 = totalSO4_[celli];

    if (totalConc_Na < 0.0)
    {
        if (totalConc_Na < -1.0 * effectiveConc_)
        {
            validTotalConc = false;
        }
        else
        {
            totalConc_Na = 0.0;
        }
    }

    if (totalConc_SO4 < 0.0)
    {
        if (totalConc_SO4 < -1.0 * effectiveConc_)
        {
            validTotalConc = false;
        }
        else
        {
            totalConc_SO4 = 0.0;
        }
    }

    if (cationTotalConc > effectiveConc_ && validTotalConc)
    {
        if (totalConc_NH3 > effectiveConc_)
        {
            equilConc = NH3_[celli];
            if (equilConc > effectiveConc_ && equilConc < totalConc_NH3)
            {
                equilConcs[indexNH3_] = equilConc;
            }
            else
            {
                validEquilConc = false;
            }

            equilConc = OH_[celli];
            if (equilConc > 1e-7)
            {
                equilConcs[indexOH_] = equilConc;
            }
            else
            {
                validEquilConc = false;
            }

            if (validEquilConc)
            {
                pConcs = -1.0*Foam::log10(equilConcs);
            }
            else
            {
                if (totalConc_NH3 > 10*cationTotalConc)
                {
                    for (i=0; i<nMetals_; i++)
                    {
                        pConcs[i] = -1.0*Foam::log10(totalConcs[i]/100);
                    }

                    pConcs[indexNH3_] = -1.0*Foam::log10(totalConc_NH3 - cationTotalConc);
                }
                else if (totalConc_NH3 > cationTotalConc)
                {
                    for (i=0; i<nMetals_; i++)
                    {
                        pConcs[i] = -1.0*Foam::log10(totalConcs[i]/10);
                    }

                    pConcs[indexNH3_] = -1.0*Foam::log10(totalConc_NH3 - cationTotalConc / 2.0);
                }
                else
                {
                    for (i=0; i<nMetals_; i++)
                    {
                        pConcs[i] = -1.0*Foam::log10(totalConcs[i]);
                    }

                    pConcs[indexNH3_] = -1.0*Foam::log10(totalConc_NH3);
                }

                conc_OH = totalConc_Na + 2*(cationTotalConc - totalConc_SO4);
                if (conc_OH > 1e-7)
                {
                    pConcs[indexOH_] = -1.0*Foam::log10(conc_OH);
                }
                else
                {
                    pConcs[indexOH_] = -1.0*Foam::log10(0.0001);
                }
            }

            for (i=0; i<nComps_; i++)
            {
                negf[i] = 0.0;
            }

            for (i=0; i<nJ; i++)
            {
                J[i] = 0.0;
            }

            totalConcs[indexNH3_] = totalConc_NH3;
            totalConcs[indexNa_] = totalConc_Na;
            totalConcs[indexSO4_] = totalConc_SO4;

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

            conc_NH3 = Foam::pow(10, -1.0*pConcs[indexNH3_]);
            conc_OH = Foam::pow(10, -1.0*pConcs[indexOH_]);

            conc_NH4 = Kb_NH3_*conc_NH3 / conc_OH;
            conc_H = kw_ / conc_OH;

            List<scalar> cationMolalConc(nTotalConc_);
            for (j=0; j<nMetals_; j++)
            {
                cationMolalConc[j] = totalConcs[j];
            }
            cationMolalConc[nMetals_] = totalConc_Na;
            cationMolalConc[nMetals_ + 1] = conc_NH4;
            cationMolalConc[nMetals_ + 2] = conc_H;

            List<scalar> anionMolalConc({totalConc_SO4, conc_OH});

            scalar I_s = activityModel_->ionic_strength(
                cationMolalConc, anionMolalConc);

            k_sp_NMC = 1.0;
            powConcs_NMC = 1.0;
            for (j=0; j<nMetals_; j++)
            {
                equilConc = Foam::pow(10, -1.0*pConcs[j]);
                metalCationConcs_[j][celli] = equilConc;

                cationConcRatio = totalConcs[j] / cationTotalConc;
                cationConcRatios_[j][celli] = cationConcRatio;

                k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);

                scalar gamma_ca = activityModel_->pair_activity_coeff(
                    j, cationMolalConc, anionMolalConc, I_s);

                powConcs_NMC *=
                    Foam::pow
                    (
                        equilConc * Foam::pow(gamma_ca, 3),
                        cationConcRatio
                    );
            }

            NH3_[celli] = conc_NH3;
            OH_[celli] = conc_OH;
            pH_[celli] = 14 - pConcs[indexOH_];

            superSat_[celli] =
                Foam::pow
                (
                    powConcs_NMC*conc_OH*conc_OH / k_sp_NMC,
                    1.0/3.0
                );
        }
        else // totalConc_NH3 <= effectiveConc_
        {
            NH3_[celli] = Foam::max(totalConc_NH3, 0.0);

            conc_OH = totalConc_Na + 2*(cationTotalConc - totalConc_SO4);

            if (conc_OH > 1e-7)
            {
                conc_NH4 = 0.0;
                conc_H = kw_ / conc_OH;

                List<scalar> cationMolalConc(nTotalConc_);
                for (j=0; j<nMetals_; j++)
                {
                    cationMolalConc[j] = totalConcs[j];
                }
                cationMolalConc[nMetals_] = totalConc_Na;
                cationMolalConc[nMetals_ + 1] = conc_NH4;
                cationMolalConc[nMetals_ + 2] = conc_H;

                List<scalar> anionMolalConc({totalConc_SO4, conc_OH});

                scalar I_s = activityModel_->ionic_strength(
                    cationMolalConc, anionMolalConc);

                k_sp_NMC = 1.0;
                powConcs_NMC = 1.0;
                for (j=0; j<nMetals_; j++)
                {
                    cationConcRatio = totalConcs[j] / cationTotalConc;
                    cationConcRatios_[j][celli] = cationConcRatio;

                    k_sp_NMC *= Foam::pow(k_sp_[j], cationConcRatio);

                    scalar gamma_ca = activityModel_->pair_activity_coeff(
                        j, cationMolalConc, anionMolalConc, I_s);

                    powConcs_NMC *=
                        Foam::pow
                        (
                            totalConcs[j] * Foam::pow(gamma_ca, 3),
                            cationConcRatio
                        );
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
            else  // conc_OH < 1e-7
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
    else // cationTotalConc <= effectiveConc_ || !validTotalConc
    {
        for (j=0; j<nMetals_; j++)
        {
            metalCationConcs_[j][celli] = 0.0;
            cationConcRatios_[j][celli] = 0.0;
        }
        NH3_[celli] = 0.0; // Don't save the equilibrium concentration

        superSat_[celli] = 1.0;

        charge_inerts = totalConc_Na - 2*totalConc_SO4;

        // This is only to estimate pH without solving full equilibrium
        conc_OH =
            (
                Foam::sqrt
                (
                    4
                   *(
                        kw_
                      + Foam::max
                        (
                            Kb_NH3_*totalConc_NH3, 0.0
                        )
                    )
                    + Foam::sqr(charge_inerts)
                )
                + charge_inerts
            ) / 2.0;

        OH_[celli] = 0.0; // Don't save the equilibrium concentration
        pH_[celli] = 14 + Foam::log10(conc_OH);

        return false;
    }

    return true;
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

    DEff = Dturb + D_[indexNa_];
    makeConcTransportEq(totalNa_);

    DEff = Dturb + D_[indexSO4_];
    makeConcTransportEq(totalSO4_);
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
    correctSpecies(totalNa_, writeSummary);
    correctSpecies(totalSO4_, writeSummary);
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
