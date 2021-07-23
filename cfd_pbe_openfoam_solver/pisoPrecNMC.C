/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Application
    pisoPrecNMC

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

#include "solutionNMC.H"
#include "populationBalance.H"
#include "odeSolver/odeSolver.H"
#include "odeSolver/UserData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    Switch CFD
    (
        piso.dict().lookupOrDefault<Switch>("CFD", true)
    );

    Switch precipitation
    (
        piso.dict().lookupOrDefault<Switch>("precipitation", false)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Foam::Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        if (CFD)
        {
            // Pressure-velocity PISO corrector
            {
                #include "UEqn.H"

                // --- PISO loop
                while (piso.correct())
                {
                    #include "pEqn.H"
                }
            }

            laminarTransport.correct();
            turbulence->correct();
        }

        /* The solution of the species and PBE comes after the piso loop,
           which is done by adopting operator-splitting approach */
        if (precipitation)
        {
            pb->transport_moments();
            solution_nmc.transport_species();

            // Precipitation is included by using the operator-splitting method
            #include "precSource.H"
        }

        runTime.write();

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Foam::Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
