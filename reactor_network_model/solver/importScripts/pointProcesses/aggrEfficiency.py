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


# Aggregation efficiency abstract class
class AggrEfficiency(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def Eff(self, superSat, L1, L2, G_L_eq):
        pass


# ---------------------- Aggregation Efficiency Models ---------------------- #
class Exponential(AggrEfficiency):

    paramList = ['A']

    def __init__(self, dict):
        self.A = init_run.strToFloat(
            init_run.read_key(dict, 'A', 'aggrEfficiency'), 'A')
        self.nu = init_run.strToFloat(
            init_run.read_key(dict, 'nu', 'aggrEfficiency'), 'nu')
        self.rhoLiq = init_run.strToFloat(
            init_run.read_key(dict, 'rhoLiq', 'aggrEfficiency'), 'rhoLiq')
        self.ti = None
        self.C0 = math.sqrt(self.rhoLiq) * (self.nu**0.25)

        super(Exponential, self).__init__()

    def update(self, epsilon):

        self.ti = math.sqrt(self.nu / epsilon)

        return self.C0 * (epsilon**0.25) / math.sqrt(10**self.A)

    def L_eq(self, L1, L2):

        return L1*L2 / math.sqrt(L1**2 + L2**2 - L1*L2)

    def Eff(self, superSat, L1, L2, Db, G_Db):

        if superSat > 1.0 and G_Db > 0.0 and L1 > 1e-9 and L2 > 1e-9:
            if L1 > L2:
                r_L = L1 / L2
            else:
                r_L = L2 / L1
            sqrt_r_L = math.sqrt(r_L**2 - 1.0)
            f_lambda = 4 * (1 + r_L - sqrt_r_L) \
                / ((1.0/3.0 + r_L - sqrt_r_L)
                   - ((r_L - sqrt_r_L)**2)*(2*r_L/3 + sqrt_r_L/3))

            tc = Db / f_lambda / G_Db

            return math.exp(-1.0 * tc / self.ti)

        return 0.0


class Rational(AggrEfficiency):

    paramList = ['A']

    def __init__(self, dict):
        self.A = init_run.strToFloat(
            init_run.read_key(dict, 'A', 'aggrEfficiency'), 'A')
        self.nu = init_run.strToFloat(
            init_run.read_key(dict, 'nu', 'aggrEfficiency'), 'nu')
        self.rhoLiq = init_run.strToFloat(
            init_run.read_key(dict, 'rhoLiq', 'aggrEfficiency'), 'rhoLiq')
        self.ti = None
        self.C0 = math.sqrt(self.rhoLiq) * (self.nu**0.25)

        super(Rational, self).__init__()

    def update(self, epsilon):

        self.ti = math.sqrt(self.nu / epsilon)

        return self.C0 * (epsilon**0.25) / math.sqrt(10**self.A)

    def L_eq(self, L1, L2):

        return L1*L2 / math.sqrt(L1**2 + L2**2 - L1*L2)

    def Eff(self, superSat, L1, L2, Db, G_Db):

        if superSat > 1.0 and G_Db > 0.0 and L1 > 1e-9 and L2 > 1e-9:
            if L1 > L2:
                r_L = L1 / L2
            else:
                r_L = L2 / L1
            sqrt_r_L = math.sqrt(r_L**2 - 1.0)
            f_lambda = 4 * (1 + r_L - sqrt_r_L) \
                / ((1.0/3.0 + r_L - sqrt_r_L)
                   - ((r_L - sqrt_r_L)**2)*(2*r_L/3 + sqrt_r_L/3))

            tc = Db / f_lambda / G_Db

            return 1.0 / (tc/self.ti + 1.0)

        return 0.0
