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
