from dataclasses import dataclass
from typing import Callable

@dataclass
class BoundaryFormula:
    value: Callable[[float, float], float]

    def __call__(self, r: float, z: float):
        return self.value(r, z)

@dataclass
class BoundaryFormulaS3(BoundaryFormula):
    beta: float
