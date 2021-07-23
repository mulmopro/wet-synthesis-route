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
