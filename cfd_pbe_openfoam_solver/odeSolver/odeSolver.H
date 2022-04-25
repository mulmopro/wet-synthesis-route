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

#ifndef odeSolver_H
#define odeSolver_H

#include "IOdictionary.H"
#include "kinematicMomentumTransportModel.H"

#include "cvode.h"
#include "nvector_serial.h"
#include "sunlinsol_lapackdense.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class solutionNMC;
class populationBalance;
struct UserData;

/*---------------------------------------------------------------------------*\
                      Class odeSolver Declaration
\*---------------------------------------------------------------------------*/

class odeSolver
:
    public IOdictionary
{

    // Private data

        //- Reference to the turbulence
        const incompressible::momentumTransportModel& turbulence_;

        //- Reference to the population balance
        const populationBalance& pb_;

        //- Reference to the solution
        const solutionNMC& solution_;

        //- problem size 
        sunindextype N_;

        //- relative tolerance
        realtype relTol_;

        //- initial step size
        realtype initialStepSize_;

        //- initial step size
        realtype maxStepSize_;

        //- initial condition
        N_Vector y0_;

        //- result
        N_Vector yout_;

        //- absolute tolerances
        N_Vector absTol_;

        //- CVODE constraint N_Vector for variables
        N_Vector constraints_;

        //- Jacobian matrix
        SUNMatrix J_;

        //- linear solver
        SUNLinearSolver LS_;

        //- pointer to the cvode solver
        void *cvode_mem_;

        //- user data (the allocated memory is freed in destructor)
        UserData *aux_data_;

        //- pointers to data of y0_ and yout_
        //- No need to be destroyed explicitly in destructor
        realtype *y0_data_, *yout_data_;

        //- Exception to terminate the integration
        //  when metal concentration is low
        struct lowMetalConcException {};

        //- Exception to prevent simulation crash
        struct invalidConcException {};

    // Private member functions

public:

    // Constructors

        //- Construct from fvMesh
        odeSolver
        (
            const fvMesh& mesh,
            const incompressible::momentumTransportModel& turbulence,
            const populationBalance& pb,
            const solutionNMC& solution
        );


    //- Destructor
    virtual ~odeSolver();


    // Member Functions

        // Public member functions

            //- Initialize and run the ODE solver 
            void solve(realtype *y, realtype t0, realtype tout);

            //- Update the user data before starting integration
            void updateUserData(label);

            //- Read ode properties dictionary
            bool read();

        // Static member functions

            //- Calculate the derivatives
            static int odeEqs(
                realtype t, N_Vector y, N_Vector ydot, void *user_data);

            //- Calculate the Jacobian
            static int Jacobian(
                realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

        // Access

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
