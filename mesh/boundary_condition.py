from dataclasses import dataclass
from typing import Callable

@dataclass
class BoundaryCondition:
    element: int
    local_border: int
    value: Callable[[float, float], float]


@dataclass
class BoundaryConditionS3(BoundaryCondition):
    beta: float
