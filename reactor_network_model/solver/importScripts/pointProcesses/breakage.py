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
