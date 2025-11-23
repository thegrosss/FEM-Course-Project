import json

from mesh.point import Point
from mesh.boundary_type import BoundaryType
from mesh.area_property import AreaProperty

from mesh.boundary_formula import BoundaryFormula, BoundaryFormulaS3


class Border:
    def __init__(self, points_indices: list[int], boundary_type: BoundaryType, formula_index: int):
        self.points_indices = points_indices
        self.boundary_type = boundary_type
        self.formula_index = formula_index

class MeshParameters:
    def __init__(self):
        self.abscissa_points_count: int = 0
        self.ordinate_points_count: int = 0
        self.control_points: list[Point] = []
        self.borders: list[Border] = []
        self.boundary_formulas: list[BoundaryFormula] = []
        self.area_properties: list[AreaProperty] = []
        self.abscissa_splits: int = 0
        self.ordinate_splits: int = 0
        self.abscissa_k: float = 0.0
        self.ordinate_k: float = 0.0
        self.refinement: int = 0

    @staticmethod
    def read_json(path: str):
        with open(path, "r") as file:
            data = json.load(file)

            def parse_point(s: str):
                s = s.replace("(", "").replace(")", "").replace(" ", "")
                r, z = map(float, s.split(","))
                return Point(r, z)

            def parse_function(s: str):
                import math

                safe_dict = {
                    'math': math,
                    'exp': math.exp,
                    'sin': math.sin,
                    'cos': math.cos,
                    'tan': math.tan,
                    'log': math.log,
                    'log10': math.log10,
                    'sqrt': math.sqrt,
                    'pi': math.pi,
                    'e': math.e
                }

                if "Ubeta" in s:
                    ubeta, beta = s.split(";")
                    ubeta = ubeta.replace("Ubeta(x,y) = ", "")
                    beta = beta.replace("beta = ", "")
                    return BoundaryFormulaS3(eval(f"lambda x, y: {ubeta}", safe_dict), float(beta))
                else:
                    s = s.replace("f(x,y) = ", "")
                    return BoundaryFormula(eval(f"lambda x, y: {s}", safe_dict))

            params = MeshParameters()
            params.abscissa_points_count = data["abscissa_points_count"]
            params.ordinate_points_count = data["ordinate_points_count"]

            params.control_points = [parse_point(p) for p in data["control_points"]]

            for ap in data["area_properties"]:
                params.area_properties.append(
                    AreaProperty(lmbda=ap["lmbda"],
                                 gamma=ap["gamma"],
                                 f=parse_function(ap["f"])))

            params.borders = [
                Border(
                    b["points_indices"],
                    BoundaryType(b["boundary_type"]),
                    b["formula_index"]
                )
                for b in data["borders"]
            ]

            params.boundary_formulas = [parse_function(f) for f in data["boundary_formulas"]]

            params.abscissa_splits = data["abscissa_splits"]
            params.ordinate_splits = data["ordinate_splits"]
            params.abscissa_k = data["abscissa_k"]
            params.ordinate_k = data["ordinate_k"]
            params.refinement = data["refinement"]

            return params