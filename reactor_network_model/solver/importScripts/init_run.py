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


import os
import math
import yaml
import numpy as np
import importScripts.exceptions as exceptions
import importScripts.activity as activity
import importScripts.chemicalEquilibria as chemicalEquilibria


def initializer(caseDir):

    # Read input files
    try:
        fileName = 'caseSetup.yml'
        filePath = os.path.join(caseDir, fileName)

        dict_case = dict()
        with open(filePath, 'r') as f:
            dict_case = yaml.safe_load(f)

        fileName = 'constantDict.yml'
        filePath = os.path.join(caseDir + '//constants', fileName)

        dict_prop = dict()
        with open(filePath, 'r') as f:
            dict_prop = yaml.safe_load(f)

    except FileNotFoundError:
        raise exceptions.FileNotFoundErr(filePath)

    ###########################################################################
    # liquid density
    liquidPropDict = read_subdict('liquidProp', dict_case)
    rhoLiq = strToFloat(
        read_key(liquidPropDict, 'rhoLiq', 'liquidProp'), 'rhoLiq')
    nu = strToFloat(read_key(liquidPropDict, 'nu', 'liquidProp'), 'nu')

    ###########################################################################
    # volumes of each compartment
    volumes = np.loadtxt('react_zone_ave.txt', skiprows=2, usecols=(3))

    # the origin and destination ID of fluxes excluding inlets
    flux_IDs = np.loadtxt(
        'react_zone_flux.txt', dtype=int, skiprows=2, usecols=(1, 3))

    # Mass flowrates excluding inlets
    massFluxes = np.loadtxt('react_zone_flux.txt', skiprows=2, usecols=(4))

    # Name of inlet
    inletOrigins = np.loadtxt(
        'react_zone_feeds.txt', dtype=str, skiprows=2, usecols=1)

    # destination ID of inlets
    inletDestinations = np.loadtxt(
        'react_zone_feeds.txt', dtype=int, skiprows=2, usecols=3)

    # Mass flow rate of the inlets
    inletMassFlowrates = np.loadtxt(
        'react_zone_feeds.txt', skiprows=2, usecols=4)

    # turbulent dissipation rate of compartments
    epsilon = np.loadtxt('react_zone_ave.txt', skiprows=2, usecols=(1))

    # turbulent kinetic energy of compartments
    kTurb = np.loadtxt('react_zone_ave.txt', skiprows=2, usecols=(2))

    # number of compartments
    nZones = np.size(volumes, axis=0)

    # volumetric flowrates
    fluxes = massFluxes / rhoLiq
    inletFlowrates = inletMassFlowrates / rhoLiq

    # outlet fluxes
    outletOrigins_list = []
    outletFlowrates_list = []

    nFluxes = np.size(flux_IDs, axis=0)
    for j in reversed(range(nFluxes)):
        if (flux_IDs[j, 0] == flux_IDs[j, 1]):
            if (fluxes[j] > 0):
                outletOrigins_list.append(int(flux_IDs[j, 0]))
                outletFlowrates_list.append(fluxes[j])

            flux_IDs = np.delete(flux_IDs, j, axis=0)
            fluxes = np.delete(fluxes, j, axis=0)

    # number of fluxes between compartments
    nFluxes = np.size(fluxes, axis=0)

    outletFlowrates = np.array(outletFlowrates_list)
    outletOrigins = np.array(outletOrigins_list)

    # correction to match sum of outlet flowrates with that of inlet flowrates
    correction = np.sum(inletFlowrates) - np.sum(outletFlowrates)

    numOutletFluxes = np.size(outletFlowrates, axis=0)
    for j in range(numOutletFluxes):
        if abs(correction/numOutletFluxes/outletFlowrates[j]) < 5e-2:
            outletFlowrates[j] += correction / numOutletFluxes
        else:
            print('ERROR: The difference between the inlet flow rates and',
                  'outlet flow rates is not negligible.')
            exit()

    ###########################################################################
    # Read atomic masses
    atomMassDict = read_subdict('atomicMass', dict_prop)
    atomicMass = list()
    atomicMass.append(strToFloat(read_key(atomMassDict, 'Ni', 'atomicMass'),
                      'Ni'))
    atomicMass.append(strToFloat(read_key(atomMassDict, 'Mn', 'atomicMass'),
                      'Mn'))
    atomicMass.append(strToFloat(read_key(atomMassDict, 'Co', 'atomicMass'),
                      'Co'))

    ###########################################################################
    # Read the start and final times for ODE integration
    input_ = read_input('startTime', dict_case)
    startTime = strToFloat(input_, 'startTime')

    input_ = read_input('finalTimes', dict_case)
    finalTimes = np.array([
        strToFloat(finalTime, 'finalTimes') for finalTime in input_])

    input_ = read_input('timeSteps', dict_case)
    timeSteps = np.array([
        strToFloat(timeStep, 'timeSteps') for timeStep in input_])

    ###########################################################################
    # Read the ode solver settings
    solverDict = read_subdict('odeSolver', dict_case)
    # solverMethod = read_key(solverDict, 'method', 'odeSolver')
    relTol = strToFloat(read_key(solverDict, 'rtol', 'odeSolver'), 'rtol')
    # solverOpts = {"method": solverMethod, "rtol": relTol}
    solverOpts = {"rtol": relTol}

    absTols = read_key(solverDict, 'atol', 'odeSolver')

    if isinstance(absTols, list):
        solverOpts["atol"] = np.array(
            [strToFloat(absTol, 'atol') for absTol in absTols])
    else:
        solverOpts["atol"] = strToFloat(absTols, 'atol')

    input_ = solverDict.pop('first_step', None)
    if input_:
        initialTimeStep = strToFloat(input_, 'first_step')
        solverOpts["h0"] = initialTimeStep

    input_ = solverDict.pop('max_step', None)
    if input_:
        maxTimeStep = strToFloat(input_, 'max_step')
        solverOpts["hmax"] = maxTimeStep

    # Read the Newton-Raphson solver settings
    nRSolverDict = read_subdict('newtonRaphson', dict_case)
    maxIter = strToInt(read_key(nRSolverDict, 'maxIter', 'newtonRaphson'),
                       'maxIter')
    tol = strToFloat(read_key(nRSolverDict, 'tol', 'newtonRaphson'), 'tol')
    nRSolverOpts = {"maxIter": maxIter, "tol": tol}

    ###########################################################################
    # Read number of nodes
    pbmDict = read_subdict('PBM', dict_case)
    nNodes = strToInt(read_key(pbmDict, 'numOfNodes', 'PBM'), 'numOfNodes')

    ###########################################################################
    # Dictionary with the input type and respective concentrations
    feedsDict = read_subdict('feeds', dict_case)

    for key in feedsDict:
        feedsDict[key]['moments'] = np.zeros(2*nNodes)

    # Read the activity model from the caseSetup file
    activityDict = read_subdict('activityModel', dict_case)
    activityModelName = read_key(activityDict, 'type', 'activityModel')
    activityModel = getattr(activity, activityModelName)()

    compNames = ['Ni', 'Mn', 'Co', 'NH3', 'Na', 'SO4']
    mixedConc = np.zeros(len(compNames))
    # flowrate-weighted average concentrations
    for i, _ in enumerate(compNames):
        for inletFlowrate, feed in zip(inletFlowrates, inletOrigins):
            mixedConc[i] += inletFlowrate * feedsDict[feed]['concentration'][i]

    # concentration of inlet species if mixed without precipitation
    mixedConc = mixedConc / (np.sum(inletFlowrates))

    cationTotalConc = np.sum(feedsDict['metals']['concentration'][0:3])
    cationConcRatios = np.zeros(3)
    aMassCrystal = 0.0
    for i, _ in enumerate(cationConcRatios):
        cationConcRatios[i] = \
            feedsDict['metals']['concentration'][i] / cationTotalConc
        aMassCrystal += \
            cationConcRatios[i] * (atomicMass[i] + 2 * (15.999 + 1.00784))

    equilibria = chemicalEquilibria.ChemicalEquilibria(
        cationConcRatios, nRSolverOpts, activityModel)

    # Read the effective concentration
    input_ = read_input('EffectiveConcentration', dict_case)
    effConc = strToFloat(input_, 'EffectiveConcentration')

    # Read crystal properties
    crystPropDict = read_subdict('crystalProp', dict_case)
    rhoCrystal = strToFloat(
        read_key(crystPropDict, 'density', 'crystalProp'), 'density')

    # Read micromixing parameter
    micromixDict = read_subdict('micromixing', dict_case)
    mixing_corr = strToFloat(
        read_key(micromixDict, 'correction_factor', 'micromixing'),
        'correction_factor')

    # print('Outlet compartments = ', outletOrigins)

    return equilibria, startTime, finalTimes, timeSteps, \
        nNodes, feedsDict, outletFlowrates, mixedConc, cationConcRatios, \
        compNames, atomicMass, aMassCrystal, rhoCrystal, solverOpts, \
        pbmDict, volumes, flux_IDs, fluxes, inletDestinations, \
        inletOrigins, inletFlowrates, outletOrigins, epsilon, kTurb, \
        nZones, nFluxes, mixing_corr, effConc, nu


# Function to check validity of the integer entered by user
def strToInt(_str, _name):
    try:
        _value = int(_str)
    except ValueError:
        print('\nError:\n"' + _name + '" should be an integer.')
        raise

    return _value


# Function to check validity of the value entered by user
def strToFloat(_str, _name):
    try:
        _value = float(_str)
    except ValueError:
        print('\nError:\nInvalid data is specified for "' + _name + '".')
        raise

    return _value


# Function to check if the required key exists
def read_key(_dict, _key, _dictName):
    try:
        _value = _dict[_key]
    except KeyError:
        raise exceptions.InputDictErr(_dictName, 'keyError', keyName=_key)

    return _value


def read_input(_key, _dict):
    try:
        _input = _dict[_key]
    except KeyError:
        raise exceptions.InputDictErr(_key, 'inputNotFound')

    return _input


def read_subdict(_subDictName, _dict):
    try:
        _subdict = _dict[_subDictName]
    except KeyError:
        raise exceptions.InputDictErr(_subDictName, 'dictNotFound')

    return _subdict


def momentLognorm(meanLog, stdLog, volFrac, nMoms):
    varLog = stdLog**2
    meanNorm = math.log((meanLog**2) / math.sqrt(varLog + meanLog**2))
    stdNorm = math.sqrt(math.log(varLog / (meanLog**2) + 1))
    varNorm = stdNorm**2

    # Number of bubbles (Moment of order 0)
    nTotal = 6.0*volFrac / math.exp(3.0*meanNorm + 9.0*varNorm/2.0) / math.pi

    moments = np.zeros(nMoms)

    for i, _ in enumerate(moments):
        moments[i] = nTotal*math.exp(i*meanNorm + 0.5*(i**2)*varNorm)

    return moments
