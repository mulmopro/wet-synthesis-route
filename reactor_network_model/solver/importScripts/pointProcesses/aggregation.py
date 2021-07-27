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
import importScripts.init_run as init_run
from abc import ABC, abstractmethod


# Aggregation abstract class
class Aggregation(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def aggregationRate(self, superSat, L1, L2):
        pass


# --------------------------- Aggregation Models ---------------------------- #

class Brownian(Aggregation):

    paramList = []
    # paramList = ['C_adj_B']

    def __init__(self, dict):
        self.C_adj_B = init_run.strToFloat(
            init_run.read_key(dict, 'C_adj_B', 'Brownian'), 'C_adj_B')
        self.kB = 1.38064852e-23  # m2 kg s-2 K-1
        self.T = init_run.strToFloat(
            init_run.read_key(dict, 'T', 'Brownian'), 'T')
        self.mu = init_run.strToFloat(
            init_run.read_key(dict, 'viscosity', 'Brownian'), 'viscosity')
        self.C0 = 2*self.kB*self.T / self.mu / 3

        super(Brownian, self).__init__()

    def aggregationRate(self, superSat, L1, L2):

        if superSat > 1.0 and (L1*L2) > 0.0:
            return (10**self.C_adj_B) * self.C0 * (L1 / L2 + L2 / L1 + 2)

        return 0.0


class Hydrodynamic(Aggregation):

    # paramList = []
    paramList = ['C_adjust']

    def __init__(self, dict):
        self.C_adjust = init_run.strToFloat(
            init_run.read_key(dict, 'C_adjust', 'Hydrodynamic'), 'C_adjust')
        self.nu = init_run.strToFloat(
            init_run.read_key(dict, 'nu', 'Hydrodynamic'), 'nu')
        self.C_aux = 2.2943 / math.sqrt(self.nu)
        self.C0 = None

        super(Hydrodynamic, self).__init__()

    def update_C0(self, epsilon):
        self.C0 = self.C_aux * math.sqrt(epsilon)

    def aggregationRate(self, superSat, L1, L2):

        if superSat > 1.0 and L1 < 2e-3 and L2 < 2e-3 and L1 > 1e-9 and L2 > 1e-9:
            return (10**self.C_adjust) * self.C0 * ((L1 + L2)**3)

        return 0.0
