from mesh.point import Point
from mesh.biquadratic_quad_element import BiquadraticQuadElement
from mesh.area_property import AreaProperty
from mesh.boundary_condition import BoundaryCondition, BoundaryConditionS3

class Mesh:
    def __init__(self,
                 points: list[Point],
                 elements: list[BiquadraticQuadElement],
                 materials: list[AreaProperty],
                 dirichlet: list[BoundaryCondition],
                 neumann: list[BoundaryCondition],
                 newton: list[BoundaryConditionS3]):

        self.points = points
        self.elements = elements
        self.materials = materials
        self.dirichlet = dirichlet
        self.neumann = neumann
        self.newton = newton