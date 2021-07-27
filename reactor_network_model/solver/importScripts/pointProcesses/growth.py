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


import importScripts.init_run as init_run
from abc import ABC, abstractmethod


# Growth abstract class
class Growth(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def growthRate(self, superSat, crystalSize, metalConc=None,
                   aMassCrystal=None, rhoCrystal=None):
        pass

    def set_ksp(self, cationConcRatios, k_sp):
        pass


# --------------------------- Growth Models --------------------------------- #

class Constant(Growth):

    paramList = ['G0']

    def __init__(self, dict):
        self.G0 = init_run.strToFloat(
            init_run.read_key(dict, 'G0', 'growth'), 'G0')
        super(Constant, self).__init__()

    def growthRate(self, superSat, crystalSize, metalConc=None,
                   aMassCrystal=None, rhoCrystal=None):

        if superSat > 1.0:
            return self.G0

        return 0.0


class PowerLaw(Growth):

    paramList = ['K', 'n']

    def __init__(self, dict):
        self.K = init_run.strToFloat(init_run.read_key(dict, 'K', 'growth'),
                                     'K')
        self.n = init_run.strToFloat(init_run.read_key(dict, 'n', 'growth'),
                                     'n')
        super(PowerLaw, self).__init__()

    def growthRate(self, superSat, crystalSize, metalConc=None,
                   aMassCrystal=None, rhoCrystal=None):

        if superSat > 1.0:
            return (10**self.K)*((superSat - 1)**self.n)

        return 0.0


class DiffusionControlled(Growth):

    # paramList = ['K', 'n']
    paramList = []

    def __init__(self, dict):
        self.k_sp_NMC = None
        self.epsilon = init_run.strToFloat(
            init_run.read_key(dict, 'epsilon', 'growth'), 'epsilon')
        self.nu = init_run.strToFloat(
            init_run.read_key(dict, 'nu', 'growth'), 'nu')
        super(DiffusionControlled, self).__init__()

    def set_ksp(self, cationConcRatios, k_sp):
        k_sp_NMC = 1
        for k_sp_i, cationConcRatio in zip(k_sp, cationConcRatios):
            k_sp_NMC *= k_sp_i**cationConcRatio
        self.k_sp_NMC = k_sp_NMC

    def growthRate(self, superSat, crystalSize, metalConc,
                   aMassCrystal, rhoCrystal):

        if superSat > 1.0 and crystalSize > 1e-20:
            metalMassConc = metalConc * 154.76  # Assuming only Nickel Sulfate
            y = (metalMassConc * 100) / (1000 + metalMassConc)
            D = -2.766e-11 * y + 6.71e-10
            Sh = 2 + 0.52*(((crystalSize**(4/3)) * (self.epsilon**(1/3))
                            / self.nu)**0.52) * ((self.nu / D)**(1/3))
            kd = Sh * D / crystalSize

            return 2*kd*((self.k_sp_NMC/4.0)**(1.0/3.0))*aMassCrystal \
                * (superSat - 1.0) / rhoCrystal

        return 0.0
