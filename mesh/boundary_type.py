from enum import Enum

class BoundaryType(Enum):
    Null = 0
    Dirichlet = 1
    Neumann = 2
    Newton = 3