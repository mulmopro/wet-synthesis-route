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


# fragment distribution abstract class
class FragmentDist(ABC):

    def __init__(self):
        super().__init__()

    @abstractmethod
    def fragment_prob(self, L, k):
        pass


# ---------------------- fragment distribution Models ----------------------- #

class symmetricFragmentation(FragmentDist):

    def __init__(self, dict):

        super(symmetricFragmentation, self).__init__()

    def fragment_prob(self, L, k):

        if L > 0.0:

            return (2**(1.0 - k/3.0)) * (L**k)

        return 0.0


class uniformFragmentation(FragmentDist):

    def __init__(self, dict):

        super(uniformFragmentation, self).__init__()

    def fragment_prob(self, L, k):

        if L > 0.0:

            return (6.0 / (k + 3.0)) * (L**k)

        return 0.0


class parabolicFragmentation(FragmentDist):

    paramList = []

    def __init__(self, dict):
        self.C = init_run.strToFloat(
            init_run.read_key(dict, 'C', 'fragmentDistribution'), 'C')

        super(parabolicFragmentation, self).__init__()

    def fragment_prob(self, L, k):

        if L > 0.0:

            return (L**k) * (
                3 * self.C / (k + 3.0) +
                (1.0 - self.C/2.0) * 18 * (6.0 - k) /
                ((k + 9.0)*(k + 6.0)*(k + 3.0)))

        return 0.0


class Laakkonen(FragmentDist):

    def __init__(self, dict):

        super(Laakkonen, self).__init__()

    def fragment_prob(self, L, k):

        if L > 0.0:

            return (L**k)*(3240.0/((k + 9.0)*(k + 12.0)*(k + 15.0)))

        return 0.0


class Erosion(FragmentDist):

    def __init__(self, dict):

        super(Erosion, self).__init__()

    def fragment_prob(self, L, k):

        if L > 0.0:

            return 1.0 + (L**3 - 1.0)**(k/3.0)

        return 0.0

