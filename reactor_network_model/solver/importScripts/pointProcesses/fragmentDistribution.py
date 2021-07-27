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
