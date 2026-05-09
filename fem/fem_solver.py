import math

from mesh.mesh import Mesh
from mesh.point import Point
from mesh.time_parameters import TimeParameters
from portrait.numerator import Numerator
from fem.basis import Basis
from fem.matrix_assembler import MatrixAssembler
from fem.los import Los
from utils import Utils


class FemSolver:
    def __init__(self, mesh: Mesh):
        Numerator.numerate_basis_functions(mesh)
        Utils.save_mesh(mesh)
        Utils.save_basis_info(mesh)

        self.mesh = mesh
        self.basis = Basis
        self.matrix_assembler = MatrixAssembler(mesh)
        self.solver = Los(10000, 1e-20)
        self.solution_time = 0.0
        self.time_layers: list[float] = []
        self.time_solutions: list[tuple[float, list[float]]] = []
        self.exact_solution = None

    def solve(self):
        matrix, vector = self.matrix_assembler.get_slae()
        matrix.print_dense("output/global_matrix")
        Utils.print_vector(vector, "global_vector")

        self.solver.compute(matrix, vector)
        Utils.save_solution(self.mesh, self.solver.solution)

    def solve_time(self, time_parameters: TimeParameters):
        initial_function = time_parameters.initial_function
        if initial_function is None:
            raise ValueError("Set exact_solution or initial_condition in input/time.json.")

        self.time_layers = time_parameters.layers
        self.exact_solution = time_parameters.exact_solution

        q_prev3 = self.build_solution_layer(initial_function, self.time_layers[0])
        q_prev2 = self.build_solution_layer(initial_function, self.time_layers[1])
        q_prev1 = self.build_solution_layer(initial_function, self.time_layers[2])

        self.time_solutions = [
            (self.time_layers[0], q_prev3.copy()),
            (self.time_layers[1], q_prev2.copy()),
            (self.time_layers[2], q_prev1.copy()),
        ]

        matrix = None
        vector = None

        for time_index in range(3, len(self.time_layers)):
            matrix, vector = self.matrix_assembler.get_hyperbolic_slae(
                self.time_layers,
                time_index,
                q_prev1,
                q_prev2,
                q_prev3,
            )

            self.solver.compute(matrix, vector)
            q = self.solver.solution.copy()

            q_prev3 = q_prev2
            q_prev2 = q_prev1
            q_prev1 = q

            self.time_solutions.append((self.time_layers[time_index], q.copy()))

        self.solver.solution = q_prev1
        self.solution_time = self.time_layers[-1]

        if matrix is not None and vector is not None:
            matrix.print_dense("output/global_matrix")
            Utils.print_vector(vector, "global_vector")

        Utils.save_solution(self.mesh, self.solver.solution)
        Utils.save_time_layers(self.time_layers)
        Utils.save_time_solutions(self.mesh, self.time_solutions)

    def get_basis_nodes(self):
        values: dict[int, Point] = {}

        for element in self.mesh.elements:
            for local_basis in range(len(element.basis_indices)):
                global_basis = element.basis_indices[local_basis]

                if global_basis not in values:
                    values[global_basis] = element.get_basis_node_position(
                        local_basis,
                        lambda idx: self.mesh.points[idx],
                    )

        return [values[i] for i in range(len(values))]

    def build_solution_layer(self, function, time: float):
        return [
            self.matrix_assembler._call_formula(function, point.r, point.z, time)
            for point in self.get_basis_nodes()
        ]

    def compare_solution_with_exact_in_nodes(self, exact_function=None, time: float | None = None):
        values: dict[int, float] = {}
        exact_function = exact_function or self.exact_solution or self.mesh.dirichlet[0].value
        time = self.solution_time if time is None else time

        for element in self.mesh.elements:
            for local_basis in range(len(element.basis_indices)):
                global_basis = element.basis_indices[local_basis]

                if global_basis not in values:
                    p = element.get_basis_node_position(local_basis, lambda idx: self.mesh.points[idx])
                    values[global_basis] = self.matrix_assembler._call_formula(exact_function, p.r, p.z, time)

        dif_square = 0.0
        exact_square = 0.0

        for i in range(len(self.solver.solution)):
            dif_square += (self.solver.solution[i] - values[i]) * (self.solver.solution[i] - values[i])
            exact_square += values[i] * values[i]

        return math.sqrt(dif_square) / math.sqrt(exact_square) if exact_square > 0.0 else math.sqrt(dif_square)

    def root_mean_square(self, points: list[Point], exact_function=None, time: float | None = None):
        func = exact_function or self.exact_solution or self.mesh.dirichlet[0].value
        time = self.solution_time if time is None else time
        dif_square = 0.0
        exact_square = 0.0

        print("Grid Exact Numeric Error")

        for p in points:
            exact = self.matrix_assembler._call_formula(func, p.r, p.z, time)
            numeric = self.value_at_point(p.r, p.z)

            dif_square += (exact - numeric) * (exact - numeric)
            exact_square += exact * exact

            print(f"{p.r:.2e} {p.z:.2e} {exact:.2e} {numeric:.2e} {abs(exact - numeric):.2e}")
        return math.sqrt(dif_square) / math.sqrt(exact_square) if exact_square > 0.0 else math.sqrt(dif_square)

    def value_at_point(self, x: float, y: float):
        point = Point(x, y)
        ielem = self.find_number_element(point)

        if ielem == -1:
            return -float('inf')

        element = self.mesh.elements[ielem]
        result = 0.0

        for i in range(9):
            global_index = element.get_global_basis_index(i)
            result += self.basis.psi(self.mesh, ielem, i, point.r, point.z) * self.solver.solution[global_index]

        return result

    def find_number_element(self, point: Point):
        for ielem in range(len(self.mesh.elements)):
            element = self.mesh.elements[ielem]
            left_bottom = self.mesh.points[element.physical_nodes_indices[0]]
            right_top = self.mesh.points[element.physical_nodes_indices[-1]]

            if left_bottom.r <= point.r <= right_top.r and left_bottom.z <= point.z <= right_top.z:
                return ielem

        return -1
