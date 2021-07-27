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


# NucleateSize abstract class
class NucleateSize(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def nucSize(self, superSat):
        pass


# --------------------------- Nucleate size Models -------------------------- #

class Fixed(NucleateSize):

    paramList = []  # ['Xc']

    def __init__(self, dict):
        self.Xc = init_run.strToFloat(
            init_run.read_key(dict, 'Xc', 'nucleateSize'), 'Xc')
        super(Fixed, self).__init__()

    def nucSize(self, superSat):

        return self.Xc
