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

        //- problem size 
        sunindextype N_;

        //- relative tolerance
        realtype relTol_;

        //- absolute tolerances
        N_Vector absTol_;

        //- initial step size
        realtype initialStepSize_;

        //- initial step size
        realtype maxStepSize_;

    // Private member functions

public:

    // Constructors

        //- Construct from fvMesh
        odeSolver
        (
            const fvMesh& mesh,
            const populationBalance& pb,
            const solutionNMC& solution
        );


    //- Destructor
    virtual ~odeSolver();


    // Member Functions

        // Public member functions

            //- Initialize and run the ODE solver 
            void solve(
                realtype *y, realtype t0, realtype tout, UserData *aux_data);

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
