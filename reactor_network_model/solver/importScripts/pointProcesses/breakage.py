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


# Breakage abstract class
class Breakage(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def breakageRate(self, superSat, L, epsilon):
        pass


# --------------------------- Breakage Models ---------------------------- #

class PowerLaw(Breakage):

    paramList = []

    def __init__(self, dict):
        self.Cbr = init_run.strToFloat(
            init_run.read_key(dict, 'Cbr', 'breakage'), 'Cbr')

        self.gamma = init_run.strToFloat(
            init_run.read_key(dict, 'gamma', 'breakage'), 'gamma')

        self.nu = init_run.strToFloat(
            init_run.read_key(dict, 'nu', 'breakage'), 'nu')

        super(PowerLaw, self).__init__()

    def breakageRate(self, superSat, L, epsilon):

        if superSat > 1.0 and L > 0.0:

            kolmogorov_length = ((self.nu**3) / epsilon)**0.25
            kolmogorov_time = (self.nu / epsilon)**0.5

            return self.Cbr*((L / kolmogorov_length)**self.gamma) \
                / kolmogorov_time

        return 0.0
