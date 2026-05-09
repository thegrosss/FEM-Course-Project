import json
from typing import Callable

from mesh.formula_parser import parse_formula


class TimeParameters:
    def __init__(
        self,
        time_begin: float,
        time_end: float,
        steps: int,
        k: float = 1.0,
        refinement: int = 0,
        grid_type: int = 0,
        exact_solution: Callable | None = None,
        initial_condition: Callable | None = None,
    ):
        self.time_begin = float(time_begin)
        self.time_end = float(time_end)
        self.steps = int(steps)
        self.k = float(k)
        self.refinement = int(refinement)
        self.grid_type = int(grid_type)
        self.exact_solution = exact_solution
        self.initial_condition = initial_condition
        self.layers = self.build_layers()

    @staticmethod
    def read_json(path: str):
        with open(path, "r") as file:
            data = json.load(file)

        interval = data["time_interval"]
        parameters = data["parameters"]

        exact_solution = data["exact_solution"]
        initial_condition = data["initial_condition"] if "initial_condition" in data else None

        return TimeParameters(
            time_begin=interval["time_begin"],
            time_end=interval["time_end"],
            steps=parameters["nt"],
            k=parameters["kt"],
            refinement=parameters["refinement"],
            grid_type=parameters["type"],
            exact_solution=parse_formula(exact_solution) if exact_solution else None,
            initial_condition=parse_formula(initial_condition) if initial_condition else None,
        )

    def build_layers(self):
        steps = self.steps * (2 ** self.refinement)
        if steps < 3:
            raise ValueError("Four-layer time scheme requires at least 3 time steps.")

        k = self.k
        if self.refinement > 0 and abs(k - 1.0) > 1e-14:
            k = k ** (1.0 / (2 ** self.refinement))

        if abs(k - 1.0) < 1e-14 or self.grid_type == 0:
            step = (self.time_end - self.time_begin) / steps
            return [self.time_begin + step * i for i in range(steps + 1)]

        step = (self.time_end - self.time_begin) * (1.0 - k) / (1.0 - k ** steps)
        layers = [self.time_begin]
        current = self.time_begin

        for _ in range(steps):
            current += step
            layers.append(current)
            step *= k

        layers[-1] = self.time_end
        return layers

    @property
    def initial_function(self):
        return self.initial_condition or self.exact_solution
