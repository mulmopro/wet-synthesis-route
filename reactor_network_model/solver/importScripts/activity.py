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
# import numpy as np
from abc import ABC, abstractmethod


# Activity abstract class
class Activity(ABC):

    def __init__(self):
        self.cationOrder = ['Ni', 'Mn', 'Co', 'Na', 'NH4', 'H']
        self.anionOrder = ['SO4', 'OH']

        self.cationDict = {
            'Ni': {'index': 0, 'Z': 2},
            'Mn': {'index': 1, 'Z': 2},
            'Co': {'index': 2, 'Z': 2},
            'Na': {'index': 3, 'Z': 1},
            'NH4': {'index': 4, 'Z': 1},
            'H': {'index': 5, 'Z': 1}}

        self.anionDict = {
            'SO4': {'index': 0, 'Z': 2},
            'OH': {'index': 1, 'Z': 1}}

        super().__init__()

    @abstractmethod
    def pair_activity(
            self, cation, anion, cationMolalConc, anionMolalConc, I_s):
        pass

    def ionic_strength(self, cationMolalConc, anionMolalConc):

        I_s = 0.0

        for i, molalConc in enumerate(cationMolalConc):
            cationName = self.cationOrder[i]
            I_s += molalConc * (self.cationDict[cationName]['Z']**2)

        for i, molalConc in enumerate(anionMolalConc):
            anionName = self.anionOrder[i]
            I_s += molalConc * (self.anionDict[anionName]['Z']**2)

        I_s *= 0.5

        return I_s


class IdealSolution(Activity):

    def __init__(self):

        super(IdealSolution, self).__init__()

    def pair_activity(
            self, cation, anion, cationMolalConc, anionMolalConc, I_s):

        return 1.0, 0.0


class Davies(Activity):

    def __init__(self):
        self.A_gamma = 0.511  # valid for T = 25 C
        self.B = 0.2  # also 0.3 is used in the literature

        super(Davies, self).__init__()

    def pair_activity(
            self, cation, anion, I_s):

        z_c = self.cationDict[cation]['Z']
        z_a = self.anionDict[anion]['Z']

        nu_c = z_a
        nu_a = z_c

        _, log_gamma_c = self.ion_activity(z_c, I_s)
        _, log_gamma_a = self.ion_activity(z_a, I_s)

        log_gamma = (nu_c * log_gamma_c + nu_a * log_gamma_a) / (nu_c + nu_a)

        return 10**log_gamma, log_gamma

    def ion_activity(self, z_i, I_s):

        sqrt_I = math.sqrt(I_s)

        log_gamma_ion = -self.A_gamma * z_i * z_i * \
            (sqrt_I / (1 + sqrt_I) - self.B * I_s)

        return 10**log_gamma_ion, log_gamma_ion


class Bromley(Activity):

    def __init__(self):
        self.A_gamma = 0.511  # valid for T = 25 C

        self.initialize()

        super(Bromley, self).__init__()

    def pair_activity(
            self, cation, anion, cationMolalConc, anionMolalConc, I_s):

        z_c = self.cationDict[cation]['Z']
        z_a = self.anionDict[anion]['Z']

        _, log_gamma_pure_ac = self.pure_pair_activity(
            cation, anion, z_c, z_a, I_s)

        f_c = 0.0
        f_a = 0.0

        sqrt_I = math.sqrt(I_s)

        coeff = self.A_gamma * sqrt_I / (1 + sqrt_I)

        for (ion, m_ion) in zip(self.anionOrder, anionMolalConc):

            if ion != anion:
                z_ion = self.anionDict[ion]['Z']

                _, log_gamma_pure = self.pure_pair_activity(
                    cation, ion, z_c, z_ion, I_s)

                y_ac = ((z_c + z_ion)**2) * m_ion / (4 * I_s)

                f_c += y_ac * (log_gamma_pure + coeff*z_c*z_ion)
            else:
                y_ac = ((z_c + z_a)**2) * m_ion / (4 * I_s)

                f_c += y_ac * (log_gamma_pure_ac + coeff*z_c*z_a)

        for (ion, m_ion) in zip(self.cationOrder, cationMolalConc):

            if ion != cation:
                z_ion = self.cationDict[ion]['Z']

                _, log_gamma_pure = self.pure_pair_activity(
                    ion, anion, z_ion, z_a, I_s)

                x_ac = ((z_ion + z_a)**2) * m_ion / (4 * I_s)

                f_a += x_ac * (log_gamma_pure + coeff*z_ion*z_a)
            else:
                x_ac = ((z_c + z_a)**2) * m_ion / (4 * I_s)

                f_a += x_ac * (log_gamma_pure_ac + coeff*z_c*z_a)

        nu_c = z_a
        nu_a = z_c

        log_gamma = (
            -coeff * (nu_c*(z_c**2) + nu_a*(z_a**2))
            + (nu_c*f_c + nu_a*f_a)) / (nu_c + nu_a)

        return 10**log_gamma, log_gamma

    def pure_pair_activity(self, cation, anion, z_c, z_a, I_s):

        i_c = self.cationDict[cation]['index']
        i_a = self.anionDict[anion]['index']

        B = self.coeff_B[i_c][i_a]

        nu_c = z_a
        nu_a = z_c

        zByZ = (nu_c*(z_c**2) + nu_a*(z_a**2)) / (nu_c + nu_a)

        sqrt_I = math.sqrt(I_s)

        log_gamma = -self.A_gamma * zByZ * sqrt_I / (1 + sqrt_I) + \
            (0.06 + 0.6*B) * zByZ * I_s / (1.0 + 1.5 * I_s / zByZ)**2 + B * I_s

        E = self.coeff_E[i_c][i_a]
        if E:
            if nu_c == 2 and nu_a == 2:
                alpha = 70
                log_gamma -= E*alpha*sqrt_I*(1 - math.exp(-alpha*sqrt_I))
                # print("Ion-association considered.")
            else:
                pass  # Not implemented because it is not needed currently

        return 10**log_gamma, log_gamma

    def initialize(self):

        paramDict = {
            'Ni': {'B': 0.054, 'delta': 0.21},
            'Mn': {'B': 0.037, 'delta': 0.21},
            'Co': {'B': 0.049, 'delta': 0.21},
            'Na': {'B': 0.0, 'delta': 0.028},
            'NH4': {'B': -0.042, 'delta': -0.02},
            'H': {'B': 0.0875, 'delta': 0.103},
            'SO4': {'B': 0.0, 'delta': -0.4},
            'OH': {'B': 0.076, 'delta': -1.0}}

        self.coeff_B = [
            [0.1056, self.estimate_B(paramDict['Ni'], paramDict['OH'])],
            [0.1226, self.estimate_B(paramDict['Mn'], paramDict['OH'])],
            [0.1244, self.estimate_B(paramDict['Co'], paramDict['OH'])],
            [-0.0204, 0.0747],
            [-0.0287, self.estimate_B(paramDict['NH4'], paramDict['OH'])],
            [0.0606, self.estimate_B(paramDict['H'], paramDict['OH'])]]

        self.coeff_E = [
            [0.00524, None],
            [0.00599, None],
            [0.00498, None],
            [None, None],
            [None, None],
            [None, None]]

    def estimate_B(self, cationParams, anionParams):

        return cationParams['B'] + anionParams['B'] + \
            cationParams['delta']*anionParams['delta']


if __name__ == "__main__":
    test_class = Bromley()
    test_class_d = Davies()

    cationMolalConc = [1.6, 0.2, 0.2, 1, 0, 0]
    anionMolalConc = [2, 1]
    I_s = test_class.ionic_strength(cationMolalConc, anionMolalConc)
    print('I_s: ', I_s)

    cationMolalConc = [0.8, 0.1, 0.1, 1, 0, 0]
    anionMolalConc = [1, 1]
    I_s = test_class.ionic_strength(cationMolalConc, anionMolalConc)
    print('I_s: ', I_s)

    ###########################################################################

    cation = 'Na'
    anion = 'OH'
    z_c = test_class.cationDict[cation]['Z']
    z_a = test_class.anionDict[anion]['Z']
    gamma_pure = test_class.pure_pair_activity(cation, anion, z_c, z_a, I_s)
    print('Bromley for pure NaOH: ', gamma_pure)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for pure NaOH: ', gamma)

    cation = 'Ni'
    anion = 'SO4'
    z_c = test_class.cationDict[cation]['Z']
    z_a = test_class.anionDict[anion]['Z']
    gamma_pure = test_class.pure_pair_activity(cation, anion, z_c, z_a, I_s)
    print('Bromley for pure NiSO4: ', gamma_pure)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for pure NiSO4: ', gamma)

    cation = 'Ni'
    anion = 'OH'
    z_c = test_class.cationDict[cation]['Z']
    z_a = test_class.anionDict[anion]['Z']
    gamma_pure = test_class.pure_pair_activity(cation, anion, z_c, z_a, I_s)
    print('Bromley for pure Ni(OH)2: ', gamma_pure)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for pure Ni(OH)2: ', gamma)

    cation = 'Ni'
    anion = 'OH'
    cationMolalConc = [5/3, 0.0, 0.0, 0, 0, 0]
    anionMolalConc = [0, 10/3]
    I_s = test_class.ionic_strength(cationMolalConc, anionMolalConc)
    print('I_s', I_s)
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for pure Ni(OH)2: ', gamma)

    ###########################################################################

    cation = 'Ni'
    anion = 'OH'
    cationMolalConc = [0.8, 0.1, 0.1, 1, 0, 0]
    anionMolalConc = [1, 1]
    I_s = test_class.ionic_strength(cationMolalConc, anionMolalConc)
    print('I_s', I_s)
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for Ni(OH)2: ', gamma)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for Ni(OH)2: ', gamma)

    cation = 'Mn'
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for Mn(OH)2: ', gamma)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for Mn(OH)2: ', gamma)

    cation = 'Co'
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for Co(OH)2: ', gamma)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for Co(OH)2: ', gamma)

    ###########################################################################

    cation = 'Ni'
    anion = 'OH'
    cationMolalConc = [0.08, 0.01, 0.01, 0.1, 0, 0]
    anionMolalConc = [0.1, 0.1]
    I_s = test_class.ionic_strength(cationMolalConc, anionMolalConc)
    print('I_s', I_s)
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for Ni(OH)2: ', gamma)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for Ni(OH)2: ', gamma)

    cation = 'Mn'
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for Mn(OH)2: ', gamma)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for Mn(OH)2: ', gamma)

    cation = 'Co'
    gamma = test_class.pair_activity(
        cation, anion, cationMolalConc, anionMolalConc, I_s)
    print('Bromley for Co(OH)2: ', gamma)
    gamma = test_class_d.pair_activity(cation, anion, I_s)
    print('Davies for Co(OH)2: ', gamma)
