'''----------------------------------------------------------------------------
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
----------------------------------------------------------------------------'''


import inspect  # argparse, importlib, inspect
import os
import sys
import traceback
import math
import time
import ray
import numpy as np
from contextlib import redirect_stdout
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from importScripts.momentCalc import adaptiveWheeler


# Entry to the main function
def main(argv):

    # Find location of the executed "pythermo.py" script
    frame = inspect.currentframe()
    filePath = inspect.getfile(frame)
    caseDir = os.path.realpath(os.path.abspath(os.path.dirname(filePath)))
    # Path to the import files
    importDir = os.path.join(caseDir, 'importScripts')

    # Start the calculations
    try:
        import importScripts.exceptions as exceptions
        from importScripts.init_run import initializer

        equilibria, startTime, finalTimes, timeSteps, nNodes, feedsDict, \
            outletFlowrates, mixedConc, cationConcRatios, \
            compNames, atomicMass, aMassCrystal, rhoCrystal, solverOpts, \
            pbmDict, num_cpus, volumes, flux_IDs, fluxes, inletDestinations, \
            inletOrigins, inletFlowrates, outletOrigins, epsilon, kTurb, \
            nZones, nFluxes, mixing_corr, effConc, nu = initializer(caseDir)

        # ray.init(local_mode=True)
        ray.init(num_cpus=num_cpus, num_gpus=0, include_dashboard=False)

        nSolve = finalTimes.size

        nEnv = len(feedsDict.keys())
        nConc = len(compNames)
        nMom = nNodes*2

        from importScripts.NiMnCoHydroxidePrec import compartmentModel

        precModel = compartmentModel(
            equilibria, mixedConc[5], nNodes, cationConcRatios, aMassCrystal,
            rhoCrystal, pbmDict, nZones, nFluxes, nEnv, nConc, feedsDict,
            mixing_corr, effConc, epsilon, kTurb, nu, num_cpus)

        resultPath = os.path.join(caseDir, 'timeResults')
        dir_0 = None
        if os.path.isdir(resultPath):
            timeDirs = [
                os.path.join(resultPath, item) for item in os.listdir(resultPath)
                if os.path.isdir(os.path.join(resultPath, item))]
            
            for timeDir in timeDirs:
                if float(os.path.basename(timeDir)) == startTime:
                    dir_0 = timeDir
                    break

        if dir_0:
            p = np.load(dir_0 + '/p.npy')
            concs = np.load(dir_0 + '/totalConc.npy')
            moments = np.load(dir_0 + '/moments.npy')

        elif startTime == 0:
            p = np.zeros((nZones, nEnv), dtype=float)
            concs = np.zeros((nZones, nConc), dtype=float)
            moments = np.zeros((nZones, nMom), dtype=float)

            # Fill the compartments with NH3
            concs[:, 3] = mixedConc[3]
            # p[:, 1] = 1.0

        else:
            print(
                "\nStart time directory \" " + os.path.basename(timeDir)
                + "\" is not found")
            exit(-1)

        begin_time = time.time()

        # print('initial environment probabilities = ', p)
        # print('initial compartments concentration = ', concs)

        if finalTimes[0] <= 1:
            writeInterval = finalTimes[0]
        elif finalTimes[0] <= 10000:
            writeInterval = finalTimes[0] / 10
        else:
            writeInterval = finalTimes[0] / 20

        # ODE solver
        tEnd, p, concs, moments, solStatus = precModel.solveODE(
            p, concs, moments, startTime, finalTimes[0], timeSteps[0],
            writeInterval, solverOpts, volumes, fluxes, flux_IDs,
            inletFlowrates, inletOrigins, inletDestinations,
            outletFlowrates, outletOrigins)

        for i in range(nSolve - 1):
            if solStatus == 0:
                # Update first step
                solverOpts["first_step"] = timeSteps[i + 1] / 1000

                # solverOptions["max_step"] = timeSteps[i + 1] / 10

                if finalTimes[i + 1] <= 1:
                    writeInterval = finalTimes[i + 1]
                elif finalTimes[i + 1] <= 10000:
                    writeInterval = finalTimes[i + 1] / 10
                else:
                    writeInterval = finalTimes[i + 1] / 20

                tEnd, p, concs, moments, solStatus = precModel.solveODE(
                    p, concs, moments, tEnd, finalTimes[i + 1],
                    timeSteps[i + 1], writeInterval, solverOpts, volumes,
                    fluxes, flux_IDs, inletFlowrates, inletOrigins,
                    inletDestinations, outletFlowrates, outletOrigins)

            else:
                break

        end_time = time.time()
        print("\n Computational Time (seconds)", end_time - begin_time)
        if solStatus == -1:
            raise exceptions.IntegrationFailedErr()
        elif solStatus == 2:
            print('The simulation has reached steady state. \
                Simulation completed')

        # checkMassBalance(
        #     nID, fluxes, fluxes_IDs, feedsDict, outletCompartments,
        #     outletFlowrates, inletOrigins, inletDestinations, inletFlowrates,
        #     solConcs[:, -1], nNodes*2, solMoments[:, -1], cationConcRatios,
        #     atomicMass, rhoCrystal, aMassCrystal)

    # If the necessary files to be imported are not found
    # in the ...\caseDir\importScripts
    except ImportError as e:
        print(str(e))
        print('\nThe above python script is not found \
              in the following address:\n' + importDir)
    # All other exceptions raised during calculation will be processed here
    except Exception as e:
        # List of handled exceptions:
        exceptionClasses = \
            ['InputDictErr', 'FileNotFoundErr', 'ValueError',
             'RealizabilityErr', 'IntegrationFailedErr']
        # If the exception is among handled exceptions,
        # print the prepared error message
        if e.__class__.__name__ in exceptionClasses:
            print(str(e))
            sys.exit(1)
        # If the exception is not foreseen,
        # print the exception message and traceback for debugging
        else:
            print(e.__class__.__name__)
            e = exceptions.create_unhandled_exception(e)
            print(str(e))
            traceback.print_exc()
            sys.exit(1)
    finally:
        print('\n')


def checkMassBalance(
        nID, fluxes, fluxes_IDs, feedsDict, outletCompartments,
        outletFlowrates, inletOrigins, inletDestinations,
        inletFlowrates, solConcs, nMoments, solMoments,
        cationConcRatios, atomicMass, rhoCrystal, aMassCrystal):

    # calculate the mass error for each compartment
    # numbers of components
    compNumber = int(np.size(solConcs)/nID)
    # number of metals
    mNumber = compNumber - 3
    # inlet and outlet mass for each compartment
    massConc_in = np.zeros((nID, mNumber))
    massConc_out = np.zeros((nID, mNumber))
    # Mass error for each compartmens
    massConcErrors = np.zeros((nID, mNumber))
    for ID in range(nID):
        massCrystal = np.zeros(3)
        m3_in = 0
        m3_out = 0
        for i, flux in enumerate(fluxes):
            From = fluxes_IDs[i, 0]
            To = fluxes_IDs[i, 1]
            if From == ID:
                for k, aMass in enumerate(atomicMass):
                    massConc_out[ID, k] += \
                        solConcs[compNumber * ID + k] * aMass * flux
                m3_out += solMoments[nMoments * ID + 3] * flux
            if To == ID:
                for k, aMass in enumerate(atomicMass):
                    massConc_in[ID, k] += solConcs[ID + k] * aMass * flux
                m3_in += solMoments[nMoments * ID + 3] * flux

        for outFlow, outComp in \
                zip(outletFlowrates, outletCompartments):
            if ID == outComp:
                for k, aMass in enumerate(atomicMass):
                    massConc_out[ID, k] += \
                        solConcs[compNumber * ID + k] * aMass * outFlow
                m3_out += solMoments[nMoments * ID + 3] * outFlow

        for inletID, inletOrigin, inletFlowrate in \
                zip(inletDestinations, inletOrigins, inletFlowrates):
            if ID == inletID:
                for k, aMass in enumerate(atomicMass):
                    massConc_in[ID, k] += \
                        feedsDict[inletOrigin]['concentration'][k] * \
                        aMass * inletFlowrate

        for i in range(mNumber):
            massCrystal[i] = \
                (m3_out - m3_in) * rhoCrystal * math.pi * \
                cationConcRatios[i] * atomicMass[i] / (6.0 * aMassCrystal)

        for i in range(mNumber):
            massConcErrors[ID, i] = \
                massConc_in[ID, i] - massConc_out[ID, i] - massCrystal[i]

    ##################################################################

    # calculate the mass error for the total system
    massConc_in_tot = np.zeros(mNumber)
    massConc_out_tot = np.zeros(mNumber)
    massConcErrors_tot = np.zeros(mNumber)
    massCrystal_tot = np.zeros(mNumber)
    m3_out_tot = 0
    for outComp, outFlow in \
            zip(outletCompartments, outletFlowrates):
        for k, aMass in enumerate(atomicMass):
            massConc_out_tot[k] += \
                solConcs[compNumber * outComp + k] \
                * aMass * outFlow
        m3_out_tot += \
            solMoments[nMoments * outComp + 3] * outFlow

    for inletOrigin, inletFlowrate in zip(inletOrigins, inletFlowrates):
        for k, aMass in enumerate(atomicMass):
            massConc_in_tot[k] += \
                feedsDict[inletOrigin]['concentration'][k] * \
                aMass * inletFlowrate

    for i, (cationConcRatio, aMass) in \
            enumerate(zip(cationConcRatios, atomicMass)):
        massCrystal_tot[i] = \
            m3_out_tot * rhoCrystal * math.pi * \
            cationConcRatio * aMass / (6.0 * aMassCrystal)

    for i in range(mNumber):
        massConcErrors_tot[i] = \
            massConc_in_tot[i] - massConc_out_tot[i] - massCrystal_tot[i]

    print('total error', massConcErrors_tot)
    print('error of compartments', massConcErrors)


# Python starts the program from here because it is the only executable code
# at indentation level 0.
# Python checks if the code is executed as a script
#  and not imported by other script
if __name__ == "__main__":
    # with open('log.txt', 'w') as f:
    #     with redirect_stdout(f):
            # If it is executed as script,
            # it calls function "main" with arguments entered by the user
    main(sys.argv[1:])
