"""Distribution classes for indel length modeling"""

import _Sailfish
import math
from typing import List

class Distribution:
    """Base class for discrete distributions"""
    
    def set_dist(self, dist):
        epsilon = 10e-6
        if abs(sum(dist)-1) > epsilon:
            raise ValueError(f"Sum must be 1, got {sum(dist)}")
        for x in dist:
            if x < 0 or x > 1:
                raise ValueError(f"Values must be in [0,1], got {x}")
        self._dist = _Sailfish.DiscreteDistribution(dist)
    
    def _get_Sailfish_dist(self) -> _Sailfish.DiscreteDistribution:
        return self._dist
    

class CustomDistribution(Distribution):
    '''
    Provide a custom discrete distribution to the model.
    '''
    def __init__(self, dist: List[float]):
        self.set_dist(dist)

class GeometricDistribution(Distribution):
    def __init__(self, p: float, truncation: int = 150):
        """
        Calculation of geoemtric moment
        inputs:
        p - p parameter of the geoemtric distribution
        truncation - (optional, by default 150) maximal value of the distribution
        """
        self.p = p
        self.truncation = truncation
        def PMF(x):
            return p*(1-p)**(x-1)
        
        def CDF(x):
            return 1-(1-p)**x
        
        norm_factor = CDF(truncation) - CDF(0)

        probabilities = [PMF(i)/norm_factor for i in range(1, truncation+1)]
        # probabilities = probabilities / norm_factor

        self.set_dist(probabilities)

    def __repr__(self) -> str:
        return f"Geometric distribution: (p={self.p}, truncation{self.truncation})"

class PoissonDistribution(Distribution):
    def __init__(self, p: float, truncation: int = 150):
        """
        Calculation of geoemtric moment
        inputs:
        p - p parameter of the geoemtric distribution
        truncation - (optional, by default 150) maximal value of the distribution
        """
        self.p = p
        self.truncation = truncation
        factorial = math.factorial
        # factorial = lambda z: reduce(operator.mul, [1, 1] if z == 0 else range(1,z+1))

        def PMF(x):
            return ((p**x)*(math.e**-p))*(1.0/factorial(x))
        def CDF(x):
            return (math.e**-p)*sum([(p**i)*(1.0/factorial(i)) for i in range(0,x+1)])

        norm_factor = CDF(truncation) - CDF(0)

        probabilities = [PMF(i)/norm_factor for i in range(1, truncation+1)]

        self.set_dist(probabilities)

    def __repr__(self) -> str:
        return f"Poisson distribution: (p={self.p}, truncation{self.truncation})"

class ZipfDistribution(Distribution):
    def __init__(self, p: float, truncation: int = 150):
        """
        Calculation of geoemtric moment
        inputs:
        p - p parameter of the geoemtric distribution
        truncation - (optional, by default 150) maximal value of the distribution
        """
        self.p = p
        self.truncation = truncation

        norm_factor = sum([(i**-p) for i in range(1,truncation+1)])
        probabilities = [(i**-p)/norm_factor for i in range(1, truncation+1)]

        self.set_dist(probabilities)
    
    def __repr__(self) -> str:
        return f"Zipf distribution: (p={self.p}, truncation{self.truncation})" 

