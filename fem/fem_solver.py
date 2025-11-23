import math

from mesh.mesh import Mesh
from mesh.point import Point
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

    def solve(self):
        matrix, vector = self.matrix_assembler.get_slae()
        matrix.print_dense("output/global_matrix")
        Utils.print_vector(vector, "global_vector")

        self.solver.compute(matrix, vector)
        Utils.print_vector(self.solver.solution, "solution")

    def compare_solution_with_exact_in_nodes(self):
        values: dict[int, float] = {}
        exact_function = self.mesh.dirichlet[0].value

        for element in self.mesh.elements:
            for local_basis in range(len(element.basis_indices)):
                global_basis = element.basis_indices[local_basis]

                if not global_basis in values:
                    p = element.get_basis_node_position(local_basis, lambda idx: self.mesh.points[idx])
                    values[global_basis] = exact_function(p.r, p.z)

        dif_square = 0.0
        exact_square = 0.0

        for i in range(len(self.solver.solution)):
            dif_square += (self.solver.solution[i] - values[i]) * (self.solver.solution[i] - values[i])
            exact_square += values[i] * values[i]

        return math.sqrt(dif_square) / math.sqrt(exact_square)

    def root_mean_square(self, points: list[Point]):
        func = self.mesh.dirichlet[0].value
        dif_square = 0.0
        exact_square = 0.0

        print("Сетка Точное Численное Вектор погрешности")

        for p in points:
            exact = func(p.r, p.z)
            numeric = self.value_at_point(p.r, p.z)

            dif_square += (exact - numeric) * (exact - numeric)
            exact_square += exact * exact

            print(f"{p.r:.2e} {p.z:.2e} {exact:.2e} {numeric:.2e} {abs(exact - numeric):.2e}")
        return math.sqrt(dif_square) / math.sqrt(exact_square)

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

            if left_bottom.r <= point.r <= right_top.r and \
                left_bottom.z <= point.z <= right_top.z:
                return ielem

        return -1