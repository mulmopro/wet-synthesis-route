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

#include "odeSolver.H"

#include "UserData.H"
#include "solutionNMC.H"
#include "populationBalance.H"
#include "growthModel.H"
#include "nucleationRateModel.H"
#include "nucleateSizeModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int Foam::odeSolver::odeEqs
(
    realtype t, N_Vector y, N_Vector ydot, void *user_data
)
{
    int i;
    label nMetals, nMoments;
    realtype superSat, growthRate, nucRate, nucSize, precRate, kv, crystalRho,
        crystalMW;
    realtype cationConcRatios[3];
    realtype *y_data, *ydot_data;
    UserData *aux_data;

    y_data = N_VGetArrayPointer(y);
    ydot_data = N_VGetArrayPointer(ydot);

    aux_data = static_cast<UserData*>(user_data);

    aux_data->solution_.update(y_data, aux_data, cationConcRatios);

    superSat = aux_data->superSat;

    // Info<<"superSat ("<< aux_data->cell_id << "): " <<superSat<<endl;

    growthRate = aux_data->growth_.rate(superSat);
    nucRate = aux_data->nucRate_.rate(superSat);
    nucSize = aux_data->nucleateSize_.size(superSat);

    nMetals = aux_data->nMetals;
    nMoments = aux_data->nMoments;

    kv = aux_data->kv;
    crystalRho = aux_data->crystalRho;
    crystalMW = aux_data->crystalMW;

    // derivative of m0
    ydot_data[nMetals] = nucRate;

    if (nMoments > 2)
    {
        for(i=1; i < nMoments; i++)
        {
            ydot_data[i + nMetals] = i * growthRate * y_data[i + nMetals - 1]
                + nucRate*pow(nucSize, i);
        }

        precRate = kv * ydot_data[3 + nMetals] * crystalRho / crystalMW;
    }
    else
    {
        ydot_data[nMetals + 1] = growthRate * y_data[nMetals]
                + nucRate*nucSize;

        precRate = 
            kv
          * (
                3.0*growthRate*pow(y_data[nMetals + 1], 2) / y_data[nMetals] +
                nucRate*pow(nucSize, 3)
            )
          * crystalRho / crystalMW;
    }

    for(i=0; i < nMetals; i++)
    {
        ydot_data[i] = -1.0 * precRate * cationConcRatios[i];
    }

    // Info<<"time "<<t<<endl;

    return (0);
}


int Foam::odeSolver::Jacobian(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
    void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

    return(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::odeSolver::odeSolver
(
    const fvMesh& mesh,
    const populationBalance& pb,
    const solutionNMC& solution
)
:
    IOdictionary
    (
        IOobject
        (
            "odeSolver",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    N_(solution.nMetals() + pb.numOfMoments()),

    relTol_(readScalar(lookup("relTol"))),

    initialStepSize_(readScalar(lookup("initialStepSize"))),

    maxStepSize_(readScalar(lookup("maxStepSize")))
{
    int i;
    ITstream is = lookup("absTol");

    absTol_ = NULL;
    absTol_ = N_VNew_Serial(N_);

    realtype *absTol_data;

    absTol_data = N_VGetArrayPointer(absTol_);

    if (is.nRemainingTokens() == is.tokenIndex() + 1)
    {
        scalar absTolScalar(readScalar(is));

        for (i = 0; i < N_; i++)
        {
            absTol_data[i] = absTolScalar;
        }
    }
    else
    {
        List<scalar> absTolList(is);

        if (absTolList.size() == N_)
        {
            for (i = 0; i < N_; i++)
            {
                absTol_data[i] = absTolList[i];
            }
        }
        else
        {
            FatalErrorInFunction
                << "Number of absolute tolerances should be " << N_
                << endl << absTolList.size() << " are specified"
                << endl << endl
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::odeSolver::~odeSolver()
{
    N_VDestroy(absTol_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::odeSolver::solve(realtype *y, realtype t0, realtype tout,
    UserData *aux_data)
{
    int i, flag;
    realtype t_reached;
    realtype *y0_data, *yout_data;

    N_Vector y0;
    N_Vector yout;
    SUNMatrix J;
    SUNLinearSolver LS;
    void *cvode_mem;

    y0 = NULL;
    yout = NULL;
    J = NULL;
    LS = NULL;
    cvode_mem = NULL;

    y0 = N_VNew_Serial(N_);

    cvode_mem = CVodeCreate(CV_BDF);

    J = SUNDenseMatrix(N_, N_);

    yout = N_VNew_Serial(N_);
    // SUNMatrix A = SUNDenseMatrix(N, N);

    LS = SUNLinSol_LapackDense(yout, J);

    try
    {
        y0_data = N_VGetArrayPointer(y0);

        for (i = 0; i < N_; i++)
        {
            y0_data[i] = y[i];
        }

        flag = CVodeSetUserData(cvode_mem, aux_data);

        // Info<<"start time: "<< t0 <<endl;
        // Info<<"initial data "<< y0_data[0] <<endl;

        flag = CVodeInit(cvode_mem, odeEqs, t0, y0);

        flag = CVodeSVtolerances(cvode_mem, relTol_, absTol_);

        flag = CVodeSetLinearSolver(cvode_mem, LS, J);

        // flag = CVodeSetJacFn(cvode_mem, Jacobian);

        // long int mxsteps;
        // flag = CVodeSetMaxNumSteps(cvode_mem, mxsteps);

        flag = CVodeSetInitStep(cvode_mem, initialStepSize_);

        flag = CVodeSetMaxStep(cvode_mem, maxStepSize_);

        // Info<<"end time: "<< tout <<endl;
        
        // The same matrix y0 is used for the output
        flag = CVode(cvode_mem, tout, yout, &t_reached, CV_NORMAL);

        yout_data = N_VGetArrayPointer(yout);

        // Info<<"return flag: "<< flag <<endl;
        // Info<<"initial data "<< yout_data[0] <<endl;

        for (i = 0; i < N_; i++)
        {
            y[i] = yout_data[i];
        }
    }
    catch(...)
    {
        N_VDestroy(y0);
        N_VDestroy(yout);
        CVodeFree(&cvode_mem);
        SUNMatDestroy(J);
        SUNLinSolFree(LS);
        throw;
    }

    N_VDestroy(y0);
    N_VDestroy(yout);
    CVodeFree(&cvode_mem);
    SUNMatDestroy(J);
    SUNLinSolFree(LS);
}


bool Foam::odeSolver::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        this->lookup("initialStepSize") >> initialStepSize_;

        this->lookup("maxStepSize") >> maxStepSize_;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
