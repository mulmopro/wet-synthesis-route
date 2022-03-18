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


import inspect  # argparse, importlib
import os
import sys
import traceback
import math
# import time
from mpi4py import MPI
from decimal import Decimal
# from contextlib import redirect_stdout

# Prevent np.linalg from using all processors
import mkl
mkl.set_num_threads(1)

# numpy should be imported after setting the number of threads
import numpy as np


# Entry to the main function
def main(argv):

    comm = MPI.COMM_WORLD
    # nproc = comm.Get_size()
    my_rank = comm.Get_rank()

    # Start the calculations
    try:
        import importScripts.exceptions as exceptions

        if my_rank == 0:
            # Find location of the executed script
            frame = inspect.currentframe()
            filePath = inspect.getfile(frame)
            caseDir = os.path.realpath(
                os.path.abspath(os.path.dirname(filePath)))
            # Path to the import files
            importDir = os.path.join(caseDir, 'importScripts')

            from importScripts.init_run import initializer

            equilibria, startTime, finalTimes, timeSteps, nNodes, feedsDict, \
                outletFlowrates, mixedConc, cationConcRatios, compNames, \
                atomicMass, aMassCrystal, rhoCrystal, solverOpts, pbmDict, \
                volumes, flux_IDs, fluxes, inletDestinations, inletOrigins, \
                inletFlowrates, outletOrigins, epsilon, kTurb, nZones, \
                nFluxes, mixing_corr, effConc, nu = initializer(caseDir)

            nEnv = len(feedsDict.keys())
            nConc = len(compNames)
            nMom = nNodes*2

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

            # print('initial environment probabilities = ', p)
            # print('initial compartments concentration = ', concs)

            from importScripts.NiMnCoHydroxidePrec import compartmentModel

            precModel = compartmentModel(
                equilibria, cationConcRatios, aMassCrystal, rhoCrystal, pbmDict,
                nNodes, nEnv, nConc, mixing_corr, effConc, epsilon, kTurb, nu,
                feedsDict)

        else:
            precModel = None
            solverOpts = None
            startTime = None
            finalTimes = None
            timeSteps = None
            mixedConc = None
            nZones = None
            nEnv = None
            nConc = None
            nMom = None

        precModel = comm.bcast(precModel, root=0)
        solverOpts = comm.bcast(solverOpts, root=0)
        nZones = comm.bcast(nZones, root=0)
        nEnv = comm.bcast(nEnv, root=0)
        nConc = comm.bcast(nConc, root=0)
        nMom = comm.bcast(nMom, root=0)

        # "Bcast" is better than "bcast" for numpy arrays,
        # but this communication is not critical
        startTime = comm.bcast(startTime, root=0)
        finalTimes = comm.bcast(finalTimes, root=0)
        timeSteps = comm.bcast(timeSteps, root=0)
        mixedConc = comm.bcast(mixedConc, root=0)

        if my_rank != 0:
            p = np.empty([nZones, nEnv])
            concs = np.empty([nZones, nConc])
            moments = np.empty([nZones, nMom])

        h0 = timeSteps / 100000
        hmax = timeSteps / 10

        # Patch the first values with those specified in the caseSetup.yml file
        h0[0] = solverOpts.get("h0", None)
        hmax[0] = solverOpts.get("hmax", None)

        if my_rank == 0:
            transport_args = (
                nFluxes, volumes, fluxes, flux_IDs, inletFlowrates,
                inletOrigins, inletDestinations, outletFlowrates,
                outletOrigins, feedsDict)
        else:
            transport_args = tuple([None]*10)

        comm.Barrier()
        begin_time = MPI.Wtime()

        solStatus = 0

        for i in range(finalTimes.size):
            if solStatus == 0: # and i == 0:
                # Update first and max steps
                # For the first interval, it can be None
                if not math.isnan(h0[i]):
                    solverOpts["h0"] = h0[i]

                # For the first interval, it can be None
                if not math.isnan(hmax[i]):
                    solverOpts["hmax"] = hmax[i]

                tEnd, p, concs, moments, solStatus = solveODE(
                    comm, precModel, p, concs, moments, startTime,
                    finalTimes[i], timeSteps[i], solverOpts, nZones,
                    mixedConc[5], *transport_args)

                startTime = tEnd

            else:
                break

        comm.Barrier()
        end_time = MPI.Wtime()

        if my_rank == 0:
            print("\n Computational Time (seconds)", end_time - begin_time)

            # if solStatus == -1:
            #     raise exceptions.IntegrationFailedErr()
            # elif solStatus == 2:
            #     print('The simulation has reached steady state. \
            #         Simulation completed')

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
        if my_rank == 0:
            print('\n')


def solveODE(
    comm, precModel, p, concs, moments, t0, finalTime, deltaT, solverOpts,
        nZones, mixedSO4Conc, nFluxes, volumes, fluxes, flux_IDs,
        inletFlowrates, inletOrigins, inletDestinations, outletFlowrates,
        outletOrigins, feedsDict):

    nproc = comm.Get_size()
    my_rank = comm.Get_rank()

    nEnv = p.shape[1]
    nConc = concs.shape[1]
    nMom = moments.shape[1]

    nEnvConc = nEnv + nConc

    nY = nEnv + nConc + nMom

    # Definitions needed for transport or writing results on the disk
    # (only master processor)
    if my_rank == 0:

        if finalTime <= 1:
            writeInterval = finalTime
        elif finalTime <= 10000:
            writeInterval =finalTime / 10
        else:
            writeInterval = finalTime / 20

        writeTime = t0 + writeInterval

        t_digits = max(0, max(
            -Decimal(str(finalTime)).normalize().as_tuple().exponent,
            -Decimal(str(writeInterval)).normalize().as_tuple().exponent))

        coeffMatrix = np.zeros([nZones, nZones])

        ''' Transport is done in two steps, each over half of the time-step
            (Strang operator-splitting) '''
        volByDeltaT = volumes / (deltaT / 2.0)

        # fluxes between compartments
        for i in range(nFluxes):
            flowrate = fluxes[i]
            From = flux_IDs[i, 0]
            To = flux_IDs[i, 1]
            coeffMatrix[From, From] += flowrate
            coeffMatrix[To, From] -= flowrate

        # outlet fluxes
        for i in range(np.size(outletFlowrates, axis=0)):
            outletFlowRate = outletFlowrates[i]
            From = outletOrigins[i]
            coeffMatrix[From, From] += outletFlowRate

        for ID in range(nZones):
            coeffMatrix[ID, ID] += volByDeltaT[ID]

        # distribute the zones over processors
        chunks = [ [] for _ in range(nproc)]
        for ID in range(nZones):
            quotient = ID // nproc
            remainder = ID % nproc

            if (quotient % 2) == 0:
                chunks[remainder].append(ID)
            else:
                # Fill reverse
                chunks[nproc - remainder - 1].append(ID)

        chunks_flattened = []
        list(map(chunks_flattened.extend, chunks))

        # sortIndex = sorted(range(nZones), key=chunks_flattened.__getitem__)
        sortIndex = np.argsort(np.array(chunks_flattened))

        chunk_nZone = np.asarray([len(chunk) for chunk in chunks], dtype=int)

        chunk_offset = np.insert(np.delete(chunk_nZone.cumsum(), -1), 0, 0)

        y = np.empty((nZones, nY))
        buffer = np.empty(nZones*nY)

        supersat_all = np.empty(nZones)

    else:
        chunks = None
        chunk_nZone = np.empty(0, dtype=int)
        chunk_offset = np.empty(0, dtype=int)
        y = None
        buffer = None
        supersat_all = None

    proc_zoneIDs = None
    proc_zoneIDs = comm.scatter(chunks, root=0)

    proc_nZone = len(proc_zoneIDs)

    t = t0

    timeTol = deltaT / 1000

    solStatus = 0

    ''' Advancing in time until the end of the interval.
        Only the master processor solve the transport. '''
    iter_no = 0
    while t < finalTime - timeTol:
        t += deltaT

        ''' A flag to perform some operations periodicly'''
        periodic_flag = ((iter_no % 100) == 1)

        if my_rank == 0:
            # print('t = ', t)

            ''' First transport over half of the time-step '''
            transport(
                p, nEnv, 'probabilities', volByDeltaT, coeffMatrix,
                inletDestinations, inletOrigins, inletFlowrates, feedsDict)
            # p[p < 0.0] = 0.0
            # p[p > 1.0] = 1.0

            transport(
                concs, nConc, 'concentration', volByDeltaT, coeffMatrix,
                inletDestinations, inletOrigins, inletFlowrates, feedsDict)

            transport(
                moments, nMom, 'moment', volByDeltaT, coeffMatrix,
                [], [], [], None)

            np.concatenate((p, concs, moments), axis=1, out=y)

            buffer = y[chunks_flattened, :].flatten()

        y_star = np.empty(proc_nZone*nY)

        # Scatter the data over processors
        comm.Scatterv(
            [buffer, chunk_nZone*nY, chunk_offset*nY, MPI.DOUBLE], y_star)

        ''' Calculate supersaturation for periodic load balance '''
        if periodic_flag:

            supersat = precModel.update_supersaturation(
                y_star, proc_nZone)

            comm.Barrier()

            # Gather the results from processors
            comm.Gatherv(supersat,
                [supersat_all, chunk_nZone, chunk_offset, MPI.DOUBLE])
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

        # Variables are updated in-place
        status = precModel.integrateSource(
            proc_zoneIDs, y_star, deltaT, solverOpts)

        comm.Barrier()

        # Gather the results from processors
        comm.Gatherv(y_star,
            [buffer, chunk_nZone*nY, chunk_offset*nY, MPI.DOUBLE])

        # solStatus = comm.reduce(status, op=MPI.MAX, root=0)

        ''' Second transport over half of the time-step '''
        if my_rank == 0:
            
            y = buffer.reshape((nZones, nY))[sortIndex, :]

            p = y[:, :nEnv]
            concs = y[:, nEnv : nEnvConc]
            moments = y[:, nEnvConc:]

            transport(
                p, nEnv, 'probabilities', volByDeltaT, coeffMatrix,
                inletDestinations, inletOrigins, inletFlowrates, feedsDict)
            # p[p < 0.0] = 0.0
            # p[p > 1.0] = 1.0

            transport(
                concs, nConc, 'concentration', volByDeltaT, coeffMatrix,
                inletDestinations, inletOrigins, inletFlowrates, feedsDict)

            transport(
                moments, nMom, 'moment', volByDeltaT, coeffMatrix,
                [], [], [], None)

            # averageSO4Conc = np.sum(
            #     volumes*(concs[:, 5]
            #     + p[:, 0]*self.feedsDict['metals']['concentration'][5])
            #     ) / np.sum(volumes)

            # if (1.0 - averageSO4Conc / mixedSO4Conc) < 1e-3:
            #     solStatus = 2

            if t > writeTime or abs(t - writeTime) < deltaT/10 \
                    or abs(t - finalTime) < timeTol:
                saveSolution(t, t_digits, p, concs, moments)
                writeTime += writeInterval

            if solStatus != 0:
                saveSolution(t, t_digits, p, concs, moments)
                writeTime += writeInterval
                break

            ''' Balance the load periodically based on the supersaturation '''
            if periodic_flag:
                print("t = {:.6f}".format(t))

                supersat_sortIndex = np.flip(
                    np.argsort(supersat_all[sortIndex]))

                for ID in range(nZones):
                    quotient = ID // nproc
                    remainder = ID % nproc

                    if (quotient % 2) == 0:
                        chunks[remainder][quotient] = supersat_sortIndex[ID]
                    else:
                        # Fill reverse
                        chunks[nproc - remainder - 1][quotient] = \
                            supersat_sortIndex[ID]

                chunks_flattened = []
                list(map(chunks_flattened.extend, chunks))

                sortIndex = np.argsort(np.array(chunks_flattened))
        
        if periodic_flag:
            proc_zoneIDs = comm.scatter(chunks, root=0)

        iter_no += 1

        comm.Barrier()

    return t, p, concs, moments, solStatus


def transport(y, nY, yName, volByDeltaT, coeffMatrix,
        inletDestinations, inletOrigins, inletFlowrates, feedsDict):

    for k in range(nY):

        constantVector = y[:, k]*volByDeltaT

        for To, inletOrigin, inletFlowrate in zip(
                inletDestinations, inletOrigins, inletFlowrates):
            inlet_values = np.array(
                feedsDict[inletOrigin][yName])
            constantVector[To] += inletFlowrate*inlet_values[k]

        y[:, k] = np.linalg.solve(coeffMatrix, constantVector)


def saveSolution(t, t_digits, solProbs, solConcs, solMoments):

    timeName = str(float(round(t, t_digits))).rstrip('0').rstrip('.')

    saveDir = "timeResults/" + timeName
    os.makedirs(saveDir, exist_ok=True)

    np.save(saveDir + '/p', solProbs)
    np.save(saveDir + '/totalConc', solConcs)
    np.save(saveDir + '/moments', solMoments)


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
