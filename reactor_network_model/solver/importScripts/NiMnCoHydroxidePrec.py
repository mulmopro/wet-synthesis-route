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


import math
import os
# import time
import ray
import numpy as np
import importScripts.exceptions as exceptions
from scipy.integrate import odeint
from importScripts.momentCalc import adaptiveWheeler
from importScripts.pointProcesses import growth as gr_module
from importScripts.pointProcesses import nucleation as nuc_module
from importScripts.pointProcesses import nucleateSize as nucSize_module
from importScripts.pointProcesses import aggregation as aggr_module
from importScripts.pointProcesses import aggrEfficiency as aggrEff_module
from importScripts.pointProcesses import breakage as break_module
from importScripts.pointProcesses import fragmentDistribution as fragDist_module
from importScripts.init_run import read_key
import copy


class compartmentModel(object):

    def __init__(
        self, equilibria, mixedSO4Conc, nNodes, cationConcRatios,
            aMassCrystal, rhoCrystal, pbmDict, nZones, nFluxes, nEnv, nConc,
            feedsDict, mixing_corr, effConc, epsilon, kTurb, nu, num_cpus):

        self.equilibria = equilibria
        self.nNodes = nNodes
        self.nMom = 2*nNodes
        self.nEnv = nEnv
        self.nConc = nConc
        self.cationConcRatios = cationConcRatios
        self.nMetals = cationConcRatios.shape[0]
        self.aMassCrystal = aMassCrystal
        self.rhoCrystal = rhoCrystal
        self.nu = nu
        self.numOfCores = num_cpus
        self.zeroSize = np.finfo(float).eps
        self.kv = math.pi / 6
        self.equilConc = np.zeros(5)
        self.mixing_corr = mixing_corr
        self.feedsDict = feedsDict
        self.mixedSO4Conc = mixedSO4Conc
        self.nZones = nZones
        self.nFluxes = nFluxes
        self.effConc = effConc
        self.epsilon = epsilon
        self.kTurb = kTurb

        # if kTurb > 0.0 and epsilon > 0.0:
        reynolds_l = kTurb / np.sqrt(epsilon*nu)
        c_phi = np.array(list(map(self.cPhi, reynolds_l)))
        self.gamma = self.mixing_corr*c_phi*epsilon / kTurb / 2.0

        self.growth = None
        self.nucleation = None
        self.nucleateSize = None
        self.aggregation = None
        self.aggrEfficiency = None
        self.breakage = None
        self.daughterFunction = None
        self.createModels(pbmDict)

    def solveODE(
        self, p0, conc0, mom0, t0, finalTime, deltaT, writeInterval,
            solverOpts, volumes, fluxes, flux_IDs, inletFlowrates,
            inletOrigins, inletDestinations, outletFlowrates, outletOrigins):

        feedsDict = self.feedsDict
        nZones = self.nZones

        numOfCores = self.numOfCores

        self_ref = ray.put(self)

        chunks = list()
        for i in range(numOfCores):
            chunks.append(list())

        for ID in range(self.nZones):
            quotient = ID // numOfCores
            remainder = ID % numOfCores

            if (quotient % 2) == 0:
                chunks[remainder].append(ID)
            else:
                # Fill reverse
                chunks[numOfCores - remainder - 1].append(ID)

        # print(chunks)

        solStatus = 0

        t = t0
        timeTol = deltaT / 1000

        writeTime = t0 + writeInterval

        p_new = copy.deepcopy(p0)
        conc_new = copy.deepcopy(conc0)
        mom_new = copy.deepcopy(mom0)

        coeffMatrix = np.zeros([nZones, nZones])

        ''' Transport is done in two steps, each over half of the time-step
            (Strang operator-splitting) '''
        volByDeltaT = volumes / (deltaT / 2.0)

        # fluxes between compartments
        for i in range(self.nFluxes):
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

        iter_no = 0
        while t < finalTime - timeTol:
            t += deltaT
            # print('t = ', t)

            ''' First transport over half of the time-step '''
            self.transport(
                p_new, conc_new, mom_new, volByDeltaT, coeffMatrix,
                inletDestinations, inletOrigins, inletFlowrates, feedsDict)

            ''' Balance the integration load on the processors periodically
                based on the supersaturation '''
            if (iter_no % 100) == 0:
                print('t = ', t)
                zone_supersat = self.update_supersaturation(p_new, conc_new)

                sorted_index = np.flip(np.argsort(zone_supersat))

                for ID in range(self.nZones):
                    quotient = ID // numOfCores
                    remainder = ID % numOfCores

                    if (quotient % 2) == 0:
                        chunks[remainder][quotient] = sorted_index[ID]
                    else:
                        # Fill reverse
                        chunks[numOfCores - remainder - 1][quotient] = \
                            sorted_index[ID]
                # print(chunks)

            ''' Integrate source terms over the time-step '''
            processes = []
            for worker in range(numOfCores):
                processes.append(
                    self.integrateSource.remote(
                        self_ref, chunks[worker], p_new, conc_new, mom_new,
                        deltaT, solverOpts, t))

            results = ray.get(processes)

            for worker, result in enumerate(results):
                y_worker = result[0]
                status_worker = result[1]
                for i, ID in enumerate(chunks[worker]):
                    sol_y = y_worker[i]
                    p_new[ID, :] = sol_y[:3]
                    conc_new[ID, :] = sol_y[3:9]
                    mom_new[ID, :] = sol_y[9:]

                    sol_status = status_worker[i]
                    if sol_status != 0:
                        solStatus = sol_status

            # time.sleep(0.1)

            ''' Second transport over half of the time-step '''
            self.transport(
                p_new, conc_new, mom_new, volByDeltaT, coeffMatrix,
                inletDestinations, inletOrigins, inletFlowrates, feedsDict)

            # averageSO4Conc = np.sum(
            #     volumes*(conc_new[:, 5]
            #     + p_new[:, 0]*self.feedsDict['metals']['concentration'][5])
            #     ) / np.sum(volumes)

            # if (1.0 - averageSO4Conc/self.mixedSO4Conc) < 1e-3:
            #     solStatus = 2

            if t > writeTime or abs(t - writeTime) < deltaT/10 \
                    or abs(t - finalTime) < timeTol:
                self.saveSolution(t, p_new, conc_new, mom_new)
                writeTime += writeInterval

            if solStatus != 0:
                self.saveSolution(t, p_new, conc_new, mom_new)
                writeTime += writeInterval
                break

            iter_no += 1

        return t, p_new, conc_new, mom_new, solStatus

    def transport(
        self, p_new, conc_new, mom_new, volByDeltaT, coeffMatrix,
            inletDestinations, inletOrigins, inletFlowrates, feedsDict):

        for k in range(self.nEnv):

            constantVector = p_new[:, k]*volByDeltaT

            for To, inletOrigin, inletFlowrate in zip(
                    inletDestinations, inletOrigins, inletFlowrates):
                inlet_values = np.array(
                    feedsDict[inletOrigin]['probabilities'])
                constantVector[To] += inletFlowrate*inlet_values[k]

            p_new[:, k] = np.linalg.solve(coeffMatrix, constantVector)

        # p_new[p_new < 0.0] = 0.0
        # p_new[p_new > 1.0] = 1.0

        for k in range(self.nConc):

            constantArray = conc_new[:, k]*volByDeltaT

            for To, inletOrigin, inletFlowrate in zip(
                    inletDestinations, inletOrigins, inletFlowrates):
                inlet_values = np.array(
                    feedsDict[inletOrigin]['concentration'])
                constantArray[To] += inletFlowrate*inlet_values[k]

            conc_new[:, k] = np.linalg.solve(coeffMatrix, constantArray)

        for k in range(self.nMom):

            constantArray = mom_new[:, k]*volByDeltaT

            mom_new[:, k] = np.linalg.solve(coeffMatrix, constantArray)

    def update_supersaturation(self, p, conc):

        nZones = self.nZones
        effConc = self.effConc
        nMetals = self.nMetals
        cationConcRatios = self.cationConcRatios
        equilibria = self.equilibria
        k_sp = equilibria.k_sp

        p4 = 1 - np.sum(p, axis=1)

        superSat = np.zeros(nZones)

        for ID in range(nZones):
            reactEnv_p = p4[ID]
            reactEnv_p = 1

            weightedConcs = conc[ID, :]

            if reactEnv_p > 0.0:
                concs = weightedConcs / reactEnv_p
                if np.all(concs[0:3] > effConc):
                    if concs[3] > effConc:
                        guess = np.zeros(5)
                        for i in range(nMetals):
                            guess[i] = concs[i]

                        guess[3] = concs[3]

                        if concs[4] > 0.001:
                            guess[4] = concs[4]
                        else:
                            guess[4] = 0.001

                        superSat[ID], _, _ = equilibria.solve(concs, guess)
                    else:
                        conc_OH = concs[4]
                        k_sp_NMC = 1
                        powConcs_NMC = 1
                        for i, (k_sp_i, cationConcRatio) in \
                                enumerate(zip(k_sp, cationConcRatios)):
                            k_sp_NMC *= (k_sp_i)**cationConcRatio
                            powConcs_NMC *= concs[i]**cationConcRatio
                        superSat[ID] = \
                            (powConcs_NMC*(conc_OH**2)/k_sp_NMC)**(1.0/3.0)
                else:
                    superSat[ID] = 0.0
            else:
                superSat[ID] = 0.0

        return superSat

    @staticmethod
    @ray.remote(num_cpus=1, num_returns=1)
    def integrateSource(
            self_ref, chunk_IDs, p, conc, mom, deltaT, solverOpts, t):

        y_chunks = []
        status_chunks = []

        for ID in chunk_IDs:
            y0 = np.concatenate((
                    p[ID, :], conc[ID, :], mom[ID, :]))

            self_ref.equilConc = np.zeros(5)

            y = odeint(
                self_ref.source, y0, [0, deltaT], args=(ID, self_ref, ),
                # rtol=None, atol=None, h0=1e-7, hmax=0.0, hmin=0.0,
                **solverOpts,
                tcrit=None, ixpr=0, mxstep=0, mxhnil=0, mxordn=12, mxords=5,
                full_output=False, printmessg=False, tfirst=False)

            y_chunks.append(y[-1, :])
            status_chunks.append(0)

        return y_chunks, status_chunks

    @staticmethod
    def source(y, t, ID, self):
        # print(t)

        equilibria = self.equilibria
        cationConcRatios = self.cationConcRatios
        effConc = self.effConc
        nNodes = self.nNodes
        aMassCrystal = self.aMassCrystal
        rhoCrystal = self.rhoCrystal
        growthModel = self.growth
        epsilon = self.epsilon[ID]
        # kTurb = self.kTurb[ID]
        # gamma = self.gamma[ID]
        gamma = 0.0
        nEnv = self.nEnv
        nConc = self.nConc

        dtdy = np.zeros(y.size)

        env_p = y[:nEnv]

        # reactEnv_p = 1 - np.sum(env_p)
        reactEnv_p = 1.0
        weightedConcs = y[nEnv: nEnv + nConc]
        equilConcs = np.copy(self.equilConc)

        if reactEnv_p > 0:
            concs = weightedConcs / reactEnv_p
            if np.all(concs[0:3] > effConc):
                if concs[3] > effConc:
                    guess = np.zeros(5)
                    for i in range(self.nMetals):
                        if equilConcs[i] > 0 \
                                and equilConcs[i] < concs[i]:
                            guess[i] = equilConcs[i]
                        else:
                            guess[i] = concs[i]
                    if equilConcs[4] > 1e-7:
                        guess[4] = equilConcs[4]
                    elif concs[4] > 0.001:
                        guess[4] = concs[4]
                    else:
                        guess[4] = 0.001
                    if equilConcs[3] > 0 and equilConcs[3] < concs[3]:
                        guess[3] = equilConcs[3]
                    else:
                        guess[3] = concs[3]
                    # computing supersaturation from previous solution
                    superSat, _, equilConcs = equilibria.solve(concs, guess)
                else:
                    conc_OH = concs[4]
                    k_sp_NMC = 1
                    powConcs_NMC = 1
                    for i, (k_sp_i, cationConcRatio) in enumerate(
                            zip(equilibria.k_sp, cationConcRatios)):
                        k_sp_NMC *= (k_sp_i)**cationConcRatio
                        powConcs_NMC *= concs[i]**cationConcRatio
                    superSat = (powConcs_NMC*(conc_OH**2)/k_sp_NMC)**(1.0/3.0)
            else:
                superSat = 0

            moments = y[nConc + nEnv:] / reactEnv_p
        else:
            superSat = 0.
            moments = np.zeros(nNodes*2)
            concs = np.zeros(nConc)

        self.equilConc = np.copy(equilConcs)

        # compute weights and nodes
        if (moments[0] > 0.0 and moments[1] > 0.0):
            weights, nodes = adaptiveWheeler(moments, nNodes, t)
        else:
            weights = np.zeros(nNodes)
            nodes = np.full(nNodes, self.zeroSize)

        # print(t, ID, superSat)

        # compute growth rate
        growthRates = np.zeros(nNodes)
        totalMetalConc = np.sum(concs[0:3])
        for i, node in enumerate(nodes):
            growthRates[i] = growthModel.growthRate(
                superSat, node, totalMetalConc, aMassCrystal, rhoCrystal)
        # compute nucleation rate
        nucRate = self.nucleation.nucleationRate(superSat)
        # compute nucleate size
        critSize = self.nucleateSize.nucSize(superSat)

        aggregationDict = self.aggregation
        aggrEffModel = self.aggrEfficiency
        weightByAggregation = np.zeros((nNodes, nNodes))
        for i, (nodei, weighti) in enumerate(zip(nodes, weights)):
            for j in range(i + 1):
                nodej = nodes[j]

                L_eq = aggrEffModel.L_eq(nodei, nodej)
                Db = aggrEffModel.update(epsilon) * L_eq
                G_Db = growthModel.growthRate(
                    superSat, Db, totalMetalConc, aMassCrystal,
                    rhoCrystal)
                aggrEff_ij = aggrEffModel.Eff(
                    superSat, nodei, nodej, Db, G_Db)
                aggrRate = 0
                for kernelName in aggregationDict.keys():
                    aggrkernel = aggregationDict[kernelName]
                    if kernelName in ['Hydrodynamic']:
                        aggrkernel.update_C0(epsilon)
                        aggrRate += aggrkernel.aggregationRate(
                            superSat, nodei, nodej) * aggrEff_ij
                    else:
                        aggrRate += aggrkernel.aggregationRate(
                            superSat, nodei, nodej) * aggrEff_ij
                wByAggr_ij = weighti * weights[j] * aggrRate
                weightByAggregation[i, j] = wByAggr_ij
                if j < i:
                    weightByAggregation[j, i] = wByAggr_ij

        breakageModel = self.breakage
        fragDistModel = self.fragmentDistribution

        dmdts = np.zeros(self.nMom)
        for k, _ in enumerate(moments):
            if k > 0:
                dmdtk_g = 0
                for i, (node, weight) in enumerate(zip(nodes, weights)):
                    dmdtk_g += growthRates[i]*weight*(node**(k-1))
                dmdts[k] = k*dmdtk_g + nucRate*(critSize**k)
            else:
                dmdts[k] = nucRate

            if k != 3:
                dmdtk_a = 0.0
                dmdtk_b = 0.0
                for i, (nodei, weighti) in enumerate(zip(nodes, weights)):
                    nodeiToPow3 = nodei**3
                    nodeiToPowK = nodei**k

                    for j, nodej in enumerate(nodes):
                        dmdtk_a += weightByAggregation[i, j] \
                            * (0.5*((nodeiToPow3 + nodej**3)**(k / 3.0))
                            - nodeiToPowK)

                    dmdtk_b += weighti \
                        * breakageModel.breakageRate(
                            superSat, nodei, epsilon) \
                        * (
                            fragDistModel.fragment_prob(nodei, k)
                            - nodeiToPowK
                        )

                dmdts[k] += dmdtk_a + dmdtk_b

            dmdts[k] *= reactEnv_p

        # precipitation rate
        precRate = self.kv*dmdts[3]*rhoCrystal / aMassCrystal

        # Calculate the time derivative of the probabilities
        r_p = np.zeros(3)
        dpdts = np.zeros(3)
        for i, p in enumerate(env_p):
            if p > 0 and p < 1:
                r = gamma*p*(1 - p)
            else:
                r = 0.0
            r_p[i] = r
            dpdts[i] = -r

        dcdts = np.zeros(nConc)

        for k, cationConcRatio in enumerate(cationConcRatios):
            dcdts[k] = -1.0*precRate*cationConcRatio \
                + r_p[0] * self.feedsDict['metals']['concentration'][k]

        dcdts[3] = r_p[1] * self.feedsDict['nh3']['concentration'][3]
        dcdts[4] = r_p[2] * self.feedsDict['naoh']['concentration'][4]
        dcdts[5] = r_p[0] * self.feedsDict['metals']['concentration'][-1]

        dtdy = np.concatenate((dpdts, dcdts, dmdts))

        return dtdy

    @staticmethod
    def cPhi(reynolds_l):

        # correlation parameters
        a = [0.4093, 0.6015, 0.5851, 0.09472, -0.3903, 0.1461, -0.01604]

        log10_re_l = math.log10(reynolds_l)

        c_phi = 0
        for i, a_i in enumerate(a):
            c_phi += a_i * (log10_re_l**i)

        return c_phi

    @staticmethod
    def saveSolution(t, solProbs, solConcs, solMoments):

        timeName = str(
            round(t, -int(min(math.floor(math.log10(abs(t))), 0))))

        saveDir = "timeResults/" + timeName
        os.makedirs(saveDir, exist_ok=True)

        np.save(saveDir + '/p', solProbs)
        np.save(saveDir + '/totalConc', solConcs)
        np.save(saveDir + '/moments', solMoments)

    def createModels(self, pbmDict):

        subDictNames = ('growth', 'nucleation', 'nucleateSize', 'aggregation',
                        'aggrEfficiency', 'breakage', 'fragmentDistribution')
        modelModules = (gr_module, nuc_module, nucSize_module, aggr_module,
                        aggrEff_module, break_module, fragDist_module)
        isList = (False, False, False, True, False, False, False)
        try:
            for subDictName, modelModule, isList_i in \
                    zip(subDictNames, modelModules, isList):

                subDict = read_key(pbmDict, subDictName, 'PBM')
                if isList_i:
                    subDict = subDict if subDict else dict()
                    modelInstances = dict()
                    for modelName in subDict.keys():
                        model = getattr(modelModule, modelName)
                        modelInstances[modelName] = model(subDict[modelName])
                    setattr(self, subDictName, modelInstances)

                else:
                    modelName = read_key(subDict, 'type', subDictName)
                    model = getattr(modelModule, modelName)
                    setattr(self, subDictName, model(subDict))

        except AttributeError:
            raise exceptions.InputDictErr('PBM', 'invalidModel',
                                          line=modelName, keyName=subDictName)
