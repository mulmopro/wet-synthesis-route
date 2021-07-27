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


# Nucleation abstract class
class Nucleation(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def nucleationRate(self, superSat):
        pass


# --------------------------- Nucleation Models ----------------------------- #

class Constant(Nucleation):

    paramList = ['J0']

    def __init__(self, dict):
        self.J0 = init_run.strToFloat(
            init_run.read_key(dict, 'J0', 'nucleation'), 'J0')
        super(Constant, self).__init__()

    def nucleationRate(self, superSat):

        if superSat > 1.0:
            return self.J0

        return 0.0


class PowerLaw(Nucleation):

    paramList = ['K', 'n']

    def __init__(self, dict):
        self.K = init_run.strToFloat(
            init_run.read_key(dict, 'K', 'nucleation'), 'K')
        self.n = init_run.strToFloat(
            init_run.read_key(dict, 'n', 'nucleation'), 'n')
        super(PowerLaw, self).__init__()

    def nucleationRate(self, superSat):

        if superSat > 1.0:
            return (10**self.K)*((superSat - 1)**self.n)

        return 0.0


class VolmerWeber(Nucleation):

    paramList = ['k', 'B']

    def __init__(self, dict):
        self.k = init_run.strToFloat(
            init_run.read_key(dict, 'k', 'nucleation'), 'k')
        self.B = init_run.strToFloat(
            init_run.read_key(dict, 'B', 'nucleation'), 'B')
        super(VolmerWeber, self).__init__()

    def nucleationRate(self, superSat):

        if superSat > 1.0:
            return (10**self.k)*math.exp(-1.0*math.exp(self.B) /
                                         ((math.log(superSat))**2))

        return 0.0


class HeteroHomogeneous(Nucleation):

    paramList = ['k1', 'B1', 'k2', 'B2']

    def __init__(self, dict):
        self.k1 = init_run.strToFloat(
            init_run.read_key(dict, 'k1', 'nucleation'), 'k1')
        self.B1 = init_run.strToFloat(
            init_run.read_key(dict, 'B1', 'nucleation'), 'B1')
        self.k2 = init_run.strToFloat(
            init_run.read_key(dict, 'k2', 'nucleation'), 'k2')
        self.B2 = init_run.strToFloat(
            init_run.read_key(dict, 'B2', 'nucleation'), 'B2')
        super(HeteroHomogeneous, self).__init__()

    def nucleationRate(self, superSat):

        if superSat > 1.0:
            return (10**self.k1)*math.exp(-1.0*math.exp(self.B1) /
                                          ((math.log(superSat))**2)) + \
                    (10**self.k2)*math.exp(-1.0*math.exp(self.B2) /
                                           ((math.log(superSat))**2))

        return 0.0
