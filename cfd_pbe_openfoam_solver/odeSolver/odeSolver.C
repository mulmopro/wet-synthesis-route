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

#include "odeSolver.H"

#include "UserData.H"
#include "solutionNMC.H"
#include "populationBalance.H"
#include "growthModel.H"
#include "nucleationRateModel.H"
#include "nucleateSizeModel.H"
#include "aggregationList.H"
#include "breakageList.H"
#include "inversionAlgorithm.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int Foam::odeSolver::odeEqs
(
    realtype t, N_Vector y, N_Vector ydot, void *user_data
)
{
    int i, j;
    label nMetals, nMoments;
    realtype superSat, growthRate, nucRate, nucSize, precRate, kv, crystalRho,
        crystalMW, effectiveConc, y_data_i;
    realtype cationConcRatios[3];
    realtype *y_data, *ydot_data;
    UserData *aux_data;

    y_data = N_VGetArrayPointer(y);
    ydot_data = N_VGetArrayPointer(ydot);

    aux_data = static_cast<UserData*>(user_data);

    nMetals = aux_data->nMetals;
    nMoments = aux_data->nMoments;

    effectiveConc = aux_data->effectiveConc;

    scalar cationTotalConc(0.0);
    for (i=0; i<nMetals; i++)
    {
        y_data_i = y_data[i];

        if (!(y_data_i > 0.0))
        {
            for(j=0; j < nMetals + nMoments; j++)
            {
                ydot_data[j] = 0.0;
            }

            return (1);
        }
        cationTotalConc += y_data_i;
    }

    if (!(cationTotalConc > effectiveConc))
    {
        throw lowMetalConcException{};
    }
    else if (cationTotalConc > 2*aux_data->totalSO4)
    {
        throw invalidConcException{};
    }

    List<scalar> moments(nMoments);

    for(i=0; i < nMoments; i++)
    {
        scalar moment_i = y_data[i + nMetals];
        if (moment_i < -1.0*SMALL)
        {
            for(j=0; j < nMetals + nMoments; j++)
            {
                ydot_data[j] = 0.0;
            }

            return (1);
        }
        moments[i] = moment_i;
    }

    aux_data->solution_.update(y_data, aux_data, cationConcRatios);

    superSat = aux_data->superSat;

    // Info<<"superSat ("<< aux_data->cell_id << "): " <<superSat<<endl;

    growthRate = aux_data->growth_.rate(superSat);
    nucRate = aux_data->nucRate_.rate(superSat);
    nucSize = aux_data->nucleateSize_.size(superSat);

    kv = aux_data->kv;
    crystalRho = aux_data->crystalRho;
    crystalMW = aux_data->crystalMW;

    // derivative of m0 due to the nucleation
    ydot_data[nMetals] = nucRate;

    if (nMoments > 2)
    {
        for(i=1; i < nMoments; i++)
        {
            ydot_data[i + nMetals] = i * growthRate * moments[i - 1]
                + nucRate*pow(nucSize, i);
        }

        precRate = kv * ydot_data[3 + nMetals] * crystalRho / crystalMW;

        if (moments[3] > 1e-14)
        {
            List<List<scalar>> nodesAndWeights(aux_data->invAlgm_.inversion(moments));

            List<scalar>& nodes = nodesAndWeights[0];
            List<scalar>& weights = nodesAndWeights[1];

            PhysChemData& physChemData = aux_data->physChemData_;
            physChemData.growthRate = growthRate;

            for(i=0; i < nMoments; i++)
            {
                if (i != 3)
                {
                    ydot_data[i + nMetals] +=
                        aux_data->aggregation_.source(nodes, weights, i, physChemData)
                      + aux_data->breakage_.source(nodes, weights, i, physChemData);
                }
            }
        }
    }
    else
    {
        ydot_data[nMetals + 1] = growthRate * moments[0]
                + nucRate*nucSize;

        precRate = 
            kv
          * (
                3.0*growthRate*pow(moments[1], 2) / moments[0] +
                nucRate*pow(nucSize, 3)
            )
          * crystalRho / crystalMW;

        List<scalar> nodes(1, moments[1] / moments[0]);
        List<scalar> weights(1, moments[0]);

        PhysChemData& physChemData = aux_data->physChemData_;
        physChemData.growthRate = growthRate;

        for(i=0; i < nMoments; i++)
        {
            ydot_data[i + nMetals] +=
                aux_data->aggregation_.source(nodes, weights, i, physChemData)
              + aux_data->breakage_.source(nodes, weights, i, physChemData);
        }
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
    solutionNMC& solution
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

    y0_ = NULL;
    y0_ = N_VNew_Serial(N_);
    y0_data_ = N_VGetArrayPointer(y0_);

    yout_ = NULL;
    yout_ = N_VNew_Serial(N_);
    yout_data_ = N_VGetArrayPointer(yout_);

    constraints_ = NULL;
    constraints_ = N_VNew_Serial(N_);

    N_VConst(1.0, constraints_);

    // realtype *constraints_data;

    // constraints_data = N_VGetArrayPointer(constraints_);

    // for (i = 0; i < solution.nMetals() ; i++)
    // {
    //     constraints_data[i] = 2.0;
    // }

    // for (i = solution.nMetals(); i < N_; i++)
    // {
    //     constraints_data[i] = 0.0;
    // }

    J_ = NULL;
    J_ = SUNDenseMatrix(N_, N_);

    LS_ = NULL;
    LS_ = SUNLinSol_LapackDense(yout_, J_);

    cvode_mem_ = NULL;
    cvode_mem_ = CVodeCreate(CV_BDF);

    // use the first 
    aux_data_ = new UserData(0, solution, pb, pb.turbulence());

    int flag;
    flag = CVodeSetUserData(cvode_mem_, aux_data_);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSetUserData failed"
            << endl << exit(FatalError);}

    // integration starts always from t0=0.0 as no source term depends on time
    flag = CVodeInit(cvode_mem_, odeEqs, 0.0, y0_);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeInit failed"
            << endl << exit(FatalError);}

    flag = CVodeSetMaxOrd(cvode_mem_, 2);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSetMaxOrd failed"
            << endl << exit(FatalError);}

    flag = CVodeSVtolerances(cvode_mem_, relTol_, absTol_);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSVtolerances failed"
            << endl << exit(FatalError);}

    flag = CVodeSetLinearSolver(cvode_mem_, LS_, J_);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSetLinearSolver failed"
            << endl << exit(FatalError);}

    flag = CVodeSetInitStep(cvode_mem_, initialStepSize_);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSetInitStep failed"
            << endl << exit(FatalError);}

    flag = CVodeSetMaxStep(cvode_mem_, maxStepSize_);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSetMaxStep failed"
            << endl << exit(FatalError);}

    flag = CVodeSetMaxNumSteps(cvode_mem_, 1000);
    if (flag != CV_SUCCESS){
        FatalErrorInFunction << "CVodeSetMaxNumSteps failed"
            << endl << exit(FatalError);}

    // flag = CVodeSetConstraints(cvode_mem_, constraints_);
    // if (flag != CV_SUCCESS){
    //     FatalErrorInFunction << "CVodeSetConstraints failed"
    //         << endl << exit(FatalError);}

    // flag = CVodeSetJacFn(cvode_mem_, Jacobian);

    // long int mxsteps;
    // flag = CVodeSetMaxNumSteps(cvode_mem_, mxsteps);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::odeSolver::~odeSolver()
{
    delete aux_data_;

    N_VDestroy(absTol_);
    N_VDestroy(constraints_);
    N_VDestroy(y0_);
    N_VDestroy(yout_);
    CVodeFree(&cvode_mem_);
    SUNMatDestroy(J_);
    SUNLinSolFree(LS_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::odeSolver::solve(realtype *y, realtype t0, realtype tout)
{
    int i, flag;
    realtype t_reached;

    try
    {
        for (i = 0; i < N_; i++)
        {
            y0_data_[i] = y[i];
        }

        flag = CVodeReInit(cvode_mem_, t0, y0_);
        // if (flag != CV_SUCCESS){
        //     Info << "Error CVodeReInit" << endl; throw;}

        // Info<<"end time: "<< tout <<endl;

        flag = CVode(cvode_mem_, tout, yout_, &t_reached, CV_NORMAL);

        // if (aux_data->cell_id == )
        // {
        //     Info << "flag: "<< flag << endl;
        // }

        // Info<<"return flag: "<< flag <<endl;

        // The same matrix y is used to return the result
        for (i = 0; i < N_; i++)
        {
            y[i] = yout_data_[i];
        }
    }
    catch(lowMetalConcException)
    {
        for (i = 0; i < N_; i++)
        {
            y[i] = yout_data_[i];
        }
    }
    catch(invalidConcException)
    {}
}


void Foam::odeSolver::updateUserData
(
    Foam::label celli,
    const Foam::solutionNMC& solution,
    const Foam::populationBalance& pb,
    const Foam::incompressible::momentumTransportModel& turbulence
)
{
    aux_data_->cell_id = celli;
    aux_data_->totalNH3 = solution.totalNH3()[celli];
    aux_data_->totalNa = solution.totalNa()[celli];
    aux_data_->totalSO4 = solution.totalSO4()[celli];
    aux_data_->concNH3 = solution.NH3()[celli];
    aux_data_->concOH = solution.OH()[celli];
    aux_data_->physChemData_.epsilon =
        turbulence.epsilon()().internalField()[celli];
}


bool Foam::odeSolver::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        this->lookup("initialStepSize") >> initialStepSize_;

        this->lookup("maxStepSize") >> maxStepSize_;

        this->lookup("relTol") >> relTol_;

        int i;
        ITstream is(this->lookup("absTol"));

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
                WarningInFunction
                    << "Updating absolute tolerances failed!" << endl
                    << "Number of absolute tolerances should be " << N_
                    << endl << absTolList.size() << " are specified"
                    << endl << endl;
            }
        }

        int flag;
        flag = CVodeSVtolerances(cvode_mem_, relTol_, absTol_);
        if (flag != CV_SUCCESS){
            FatalErrorInFunction << "CVodeSVtolerances failed"
                << endl << exit(FatalError);}

        flag = CVodeSetInitStep(cvode_mem_, initialStepSize_);
        if (flag != CV_SUCCESS){
            FatalErrorInFunction << "CVodeSetInitStep failed"
                << endl << exit(FatalError);}

        flag = CVodeSetMaxStep(cvode_mem_, maxStepSize_);
        if (flag != CV_SUCCESS){
            FatalErrorInFunction << "CVodeSetMaxStep failed"
                << endl << exit(FatalError);}

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
