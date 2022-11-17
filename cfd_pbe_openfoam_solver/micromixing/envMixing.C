/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "envMixing.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::envMixing::gamma_calc() const
{
    tmp<DimensionedField<scalar, volMesh>> tx
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "gamma_mixing",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            false
        )
    );

    DimensionedField<scalar, volMesh>& x = tx.ref();

    tmp<volScalarField> k(turbulence_.k());
    tmp<volScalarField> epsilon(turbulence_.epsilon());
    tmp<volScalarField> nu(turbulence_.nu());

    forAll(x, celli)
    {
        scalar k_i = k().internalField()[celli];
        scalar epsilon_i = epsilon().internalField()[celli];

        if (k_i > 0.0 and epsilon_i > 0.0)
        {
            scalar reynolds_l =
                k_i / Foam::sqrt(epsilon_i*nu().internalField()[celli]);

            scalar log10_re_l = Foam::log10(reynolds_l);

            scalar c_phi = 0.0;

            if (reynolds_l > 0.2)
            {
                if (reynolds_l < 12853)
                {
                    for (int j=0; j<a_cphi.size(); j++)
                    {
                        c_phi += a_cphi[j] * Foam::pow(log10_re_l, j);
                    }
                }
                else
                {
                    c_phi = 2.0;
                }
            }
            else
            {
                c_phi = 0.0;
            }

            x[celli] = correctionCoeff_ * c_phi * epsilon_i / k_i / 2.0;
        }
    }

    return tx;
}


void Foam::envMixing::p_flux()
{
    tmp<DimensionedField<scalar, volMesh>> gamma = gamma_calc();
    const DimensionedField<scalar, volMesh>& gamma_ref = gamma();

    forAll(gamma_ref, celli)
    {
        scalar react_p = 1.0;

        label count_p = 0;

        List<scalar> p(nPureEnvs_);

        for (int i=0; i<nPureEnvs_; i++)
        {
            scalar p_i = pureEnvP_[i][celli];

            if (p_i < 1e-20)
            {
                p[i] = 0.0;
            }
            else if (p_i > 0.9999)
            {
                count_p = 1;
                react_p = 0.0;

                for (int j=0; j<nPureEnvs_; j++)
                {
                    p[j] = 0.0;
                    pFluxes_[j][celli] = 0.0;
                    pFluxes_sp_[j][celli] = 0.0;
                }
                p[i] = 1.0;

                break;
            }
            else
            {
                p[i] = p_i;
                react_p -= p_i;
                count_p++;
            }
        }

        scalar gamma_celli = gamma_ref[celli];
        if (count_p > 1)
        {
            for (int i=0; i<nPureEnvs_; i++)
            {
                scalar p_i = p[i];

                if (p_i < 1e-20)
                {
                    pFluxes_[i][celli] = 0.0;
                    pFluxes_sp_[i][celli] = 0.0;
                }
                else /* && p_i < 0.9999 */
                {
                    pFluxes_[i][celli] = gamma_celli * p_i * (1 - p_i);
                    pFluxes_sp_[i][celli] = gamma_celli * (1 - p_i);
                }
            }
        }
        else if (count_p == 1 && react_p > 1e-20)
        {
            for (int i=0; i<nPureEnvs_; i++)
            {
                scalar p_i = p[i];

                if (p_i < 1e-20)
                {
                    pFluxes_[i][celli] = 0.0;
                    pFluxes_sp_[i][celli] = 0.0;
                }
                else /* && p_i < 0.9999 */
                {
                    pFluxes_[i][celli] = gamma_celli * p_i * (1 - p_i);
                    pFluxes_sp_[i][celli] = gamma_celli * (1 - p_i);
                }
            }
        }
        else
        {
            for (int i=0; i<nPureEnvs_; i++)
            {
                pFluxes_[i][celli] = 0.0;
                pFluxes_sp_[i][celli] = 0.0;
            }
        }

        // reactingEnvP_[celli] = react_p;

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::envMixing::envMixing
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
            "micromixing",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phi_(phi),

    turbulence_(turbulence),

    nPureEnvs_(readInt(lookup("numPureEnv"))),

    correctionCoeff_(readScalar(lookup("correctionCoeff"))),

    turbSc_
    (
        dimensionedScalar
        (
            "turbSc_", dimless, lookup("turbulentSchmidt")
        )
    ),

    reactingEnvP_
    (
        IOobject
        (
            "p_reacting",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    )
{
    if (nPureEnvs_ != 3)
    {
        FatalErrorInFunction
            << "The implementation of the micromixing model is compatible with"
            << endl << "a reactor of three separate inlets (pure environments)"
            << endl << endl
            << exit(FatalError);
    }

    List<scalar> envConcList(lookup("componentFeedConc"));

    for (label count=0; count < envConcList.size(); count++)
    {
        envConcs_.append
        (
            new dimensionedScalar
            (
                "envConc_" + Foam::name(count),
                dimConc,
                envConcList[count]
            )
        );
    }

    for (label count=0; count < nPureEnvs_; count++)
    {
        Info<< "Reading field p" << Foam::name(count + 1) << "\n" << endl;
        pureEnvP_.append
        (
            new volScalarField
            (
                IOobject
                (
                    "p_env" + Foam::name(count + 1),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        pFluxes_.append
        (
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "r_" + Foam::name(count + 1),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless/dimTime, 0.0)
            )
        );

        pFluxes_sp_.append
        (
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "r_" + Foam::name(count + 1) + "_sp",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::envMixing::~envMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::envMixing::solve()
{
    volScalarField Dturb = turbulence_.nut() / turbSc_;

    p_flux();

    // Loop over environments
    forAll(pureEnvP_, i)
    {
        // Reference to the environment for which transport eq. will be
        // defined
        volScalarField& p = pureEnvP_[i];

        tmp<fvScalarMatrix> tpEqn =
        (
            fvm::ddt(p)
          + fvm::div(phi_, p)
          - fvm::laplacian(Dturb, p)
          + fvm::Sp(pFluxes_sp_[i], p)
        );

        tpEqn.ref().relax();

        tpEqn.ref().solve(p.name());

        // p.correctBoundaryConditions();

        Info<< p.name() << " = "
            << p.weightedAverage(mesh_.V()).value()
            << "  Min(" << p.name() << ") = " << min(p).value()
            << "  Max(" << p.name() << ") = " << max(p).value()
            << endl;

        p.max(0);
        p.correctBoundaryConditions();
    }
}


void Foam::envMixing::update_react_env()
{
    forAll(reactingEnvP_, celli)
    {
        scalar react_p = 1.0;

        for (int i=0; i<nPureEnvs_; i++)
        {
            scalar p_i = pureEnvP_[i][celli];

            if (p_i < 1e-20)
            {}
            else if (p_i > 0.9999)
            {
                react_p = 0.0;
                break;
            }
            else
            {
                react_p -= p_i;
            }
        }

        reactingEnvP_[celli] = react_p;

    }

    reactingEnvP_.correctBoundaryConditions();
}


bool Foam::envMixing::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
