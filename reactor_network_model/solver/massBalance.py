import inspect
import os
import numpy as np
from importScripts.init_run import initializer
import math


def checkConcBalance(
        solP, solConcs, solMoms, nZones, nConc, nFluxes, fluxes, flux_IDs,
        outletFlowrates, outletOrigins, inletDestinations, inletOrigins,
        inletFlowrates, feedsDict, volumes, cationConcRatios, aMassCrystal,
        rhoCrystal):
    ''' This function should be used only when the steady-state
        solution is obtained'''

    nMetals = cationConcRatios.shape[0]

    zoneConcs_crystal = np.zeros((nZones, nConc))

    for ID in range(nZones):

        moment3 = solMoms[ID, 3]

        conc_crystal = (moment3 * rhoCrystal / aMassCrystal) * (math.pi / 6.0)

        for k in range(nMetals):
            zoneConcs_crystal[ID, k] += conc_crystal * cationConcRatios[k]

    zoneConc_total = solConcs + zoneConcs_crystal

    # Concentration error of compartments
    concErrors = np.zeros((nZones, nConc))

    for i in range(nFluxes):
        flowrate = fluxes[i]

        From = flux_IDs[i, 0]
        To = flux_IDs[i, 1]

        concFlux = flowrate * zoneConc_total[From, :]

        concErrors[From, :] -= concFlux
        concErrors[To, :] += concFlux

    for i in range(np.size(outletFlowrates, axis=0)):
        outletFlowRate = outletFlowrates[i]
        From = outletOrigins[i]
        concErrors[From, :] -= outletFlowRate * zoneConc_total[From, :]

    for To, inletOrigin, inletFlowrate in zip(
            inletDestinations, inletOrigins, inletFlowrates):
        inlet_values = np.array(
            feedsDict[inletOrigin]['concentration'])
        concErrors[To, :] += inletFlowrate * inlet_values

    for ID in range(nZones):
        concErrors[ID, :] /= volumes[ID]

    ##################################################################

    concError_reactor = np.zeros(nConc)

    for i in range(np.size(outletFlowrates, axis=0)):
        outletFlowRate = outletFlowrates[i]
        From = outletOrigins[i]
        concError_reactor -= outletFlowRate * zoneConc_total[From, :]

    for To, inletOrigin, inletFlowrate in zip(
            inletDestinations, inletOrigins, inletFlowrates):
        inlet_values = np.array(
            feedsDict[inletOrigin]['concentration'])
        concError_reactor += inletFlowrate * inlet_values

    concError_reactor /= np.sum(volumes)

    print('Concentration error of compartments: ', concErrors)
    print('Concentration error of reactor: ', concError_reactor)


if __name__ == "__main__":

    frame = inspect.currentframe()
    filePath = inspect.getfile(frame)
    caseDir = os.path.realpath(os.path.abspath(os.path.dirname(filePath)))
    # Path to the import files
    importDir = os.path.join(caseDir, 'importScripts')
    equilibria, startTime, finalTimes, timeSteps, \
        nNodes, feedsDict, outletFlowrates, mixedConc, cationConcRatios, \
        compNames, atomicMass, aMassCrystal, rhoCrystal, solverOpts, \
        pbmDict, num_cpus, volumes, flux_IDs, fluxes, inletDestinations, \
        inletOrigins, inletFlowrates, outletOrigins, epsilon, kTurb, \
        nZones, nFluxes, mixing_corr, effConc, nu = initializer(caseDir)

    nConc = len(compNames)

    resultPath = os.path.join(caseDir, 'timeResults')

    timeDirs = [
        os.path.join(resultPath, item) for item in os.listdir(resultPath)
        if os.path.isdir(os.path.join(resultPath, item))]

    t_array = np.array([
        float(os.path.basename(timeDir))
        for timeDir in timeDirs])

    sort_index = np.argsort(t_array)

    lastDir = timeDirs[sort_index[-1]]

    solP = np.load(lastDir + '/p.npy')
    solConcs = np.load(lastDir + '/totalConc.npy')
    solMoms = np.load(lastDir + '/moments.npy')

    checkConcBalance(
        solP, solConcs, solMoms, nZones, nConc, nFluxes, fluxes, flux_IDs,
        outletFlowrates, outletOrigins, inletDestinations, inletOrigins,
        inletFlowrates, feedsDict, volumes, cationConcRatios, aMassCrystal,
        rhoCrystal)
