from mesh.mesh_parameters import MeshParameters, Border
from mesh.point import Point
from mesh.biquadratic_quad_element import BiquadraticQuadElement
from mesh.boundary_condition import BoundaryCondition, BoundaryConditionS3
from mesh.boundary_type import BoundaryType
from mesh.edge import Edge
from mesh.mesh import Mesh
from mesh.boundary_formula import BoundaryFormulaS3

class MeshBuilder:
    def __init__(self, mesh_parameters: MeshParameters):
        self.mesh_parameters = mesh_parameters

        self.prepare_refinement()

        total_nr = self.mesh_parameters.abscissa_splits
        total_nz = self.mesh_parameters.ordinate_splits

        self.points = [Point(0.0, 0.0) for _ in range((total_nr + 1) * (total_nz + 1))]
        self.elements = [None] * (total_nr * total_nz)
        self.ir = [0] * self.mesh_parameters.abscissa_points_count
        self.iz = [0] * self.mesh_parameters.ordinate_points_count
        self.dirichlet: list[BoundaryType] = []
        self.neumann: list[BoundaryType] = []
        self.newton: list[BoundaryType] = []

    def prepare_refinement(self):
        if self.mesh_parameters.refinement == 0:
            return

        self.mesh_parameters.abscissa_splits *= 2 ** self.mesh_parameters.refinement
        sign = 1 if self.mesh_parameters.abscissa_k >= 0 else -1
        self.mesh_parameters.abscissa_k = sign * (
                abs(self.mesh_parameters.abscissa_k) ** (1.0 / (2 ** self.mesh_parameters.refinement))
        )

        self.mesh_parameters.ordinate_splits *= 2 ** self.mesh_parameters.refinement
        sign = 1 if self.mesh_parameters.ordinate_k >= 0 else -1
        self.mesh_parameters.ordinate_k = sign * (
                abs(self.mesh_parameters.ordinate_k) ** (1.0 / (2 ** self.mesh_parameters.refinement))
        )

    def create_points(self):
        total_nr = self.mesh_parameters.abscissa_splits + 1  # Nodes count on abscissa with splits
        total_nz = self.mesh_parameters.ordinate_splits + 1  # Nodes count on ordinate with splits

        primary_nr = self.mesh_parameters.abscissa_points_count  # Primary abscissa points count
        primary_nz = self.mesh_parameters.ordinate_points_count  # Primary ordinate points count

        # Forming nodes on main horizontal lines
        ordinate_splits = 0

        for i in range(primary_nz):
            a = self.mesh_parameters.control_points[i * primary_nr]  # Start
            b = self.mesh_parameters.control_points[i * primary_nr + 1]  # End

            if self.mesh_parameters.abscissa_k < 0:
                k = -1.0 / self.mesh_parameters.abscissa_k
            else:
                k = self.mesh_parameters.abscissa_k

            if abs(k - 1.0) < 1e-14:
                h = (b.r - a.r) / self.mesh_parameters.abscissa_splits
            else:
                h = (b.r - a.r) * (1.0 - k) / (1.0 - k ** self.mesh_parameters.abscissa_splits)

            self.points[total_nr * ordinate_splits] = a

            for l in range(1, self.mesh_parameters.abscissa_splits + 1):
                x = self.points[total_nr * ordinate_splits + l - 1].r + h * (k ** (l - 1))
                t = (x - a.r) / (b.r - a.r)
                y = a.z + t * (b.z - a.z)

                self.points[total_nr * ordinate_splits + l] = Point(x, y)

            ordinate_splits += self.mesh_parameters.ordinate_splits

        # Forming nodes on main vertical lines
        abscissa_splits = 0

        for i in range(primary_nr):
            a = self.mesh_parameters.control_points[i]
            b = self.mesh_parameters.control_points[primary_nr + i]

            if self.mesh_parameters.ordinate_k < 0:
                k = -1.0 / self.mesh_parameters.ordinate_k
            else:
                k = self.mesh_parameters.ordinate_k

            if abs(k - 1.0) < 1e-14:
                h = (b.z - a.z) / self.mesh_parameters.ordinate_splits
            else:
                h = (b.z - a.z) * (1.0 - k) / (1.0 - k ** self.mesh_parameters.ordinate_splits)

            self.points[abscissa_splits] = a

            for l in range(1, self.mesh_parameters.ordinate_splits + 1):
                y = self.points[abscissa_splits + (l - 1) * total_nr].z + h * (k ** (l - 1))
                t = (y - a.z) / (b.z - a.z)
                x = a.r + t * (b.r - a.r)

                self.points[abscissa_splits + l * total_nr] = Point(x, y)

            abscissa_splits += self.mesh_parameters.abscissa_splits

        # Form inner nodes
        for i in range(1, total_nz):

            if self.mesh_parameters.abscissa_k < 0:
                k = -1.0 / self.mesh_parameters.abscissa_k
            else:
                k = self.mesh_parameters.abscissa_k

            a = self.points[i * total_nr]
            b = self.points[i * total_nr + self.mesh_parameters.abscissa_splits]

            if abs(k - 1.0) < 1e-14:
                h = (b.r - a.r) / self.mesh_parameters.abscissa_splits
            else:
                h = (b.r - a.r) * (1.0 - k) / (1.0 - k ** self.mesh_parameters.abscissa_splits)

            for l in range(1, self.mesh_parameters.abscissa_splits):
                x = self.points[i * total_nr + l - 1].r + h * (k ** (l - 1))
                t = (x - a.r) / (b.r - a.r)
                y = a.z + t * (b.z - a.z)

                self.points[i * total_nr + l] = Point(x, y)

        self.ir[1] = self.mesh_parameters.abscissa_splits
        self.iz[1] = self.mesh_parameters.ordinate_splits

    def create_elements(self):
        nx = self.mesh_parameters.abscissa_splits + 1
        ny = self.mesh_parameters.ordinate_splits + 1

        ielem = 0
        for i in range(ny - 1):
            for j in range(nx - 1):
                nodes = [
                    i * nx + j,  # левый нижний
                    i * nx + j + 1,  # правый нижний
                    (i + 1) * nx + j,  # левый верхний
                    (i + 1) * nx + j + 1  # правый верхний
                ]
                self.elements[ielem] = BiquadraticQuadElement(nodes, 0)
                ielem += 1

    def create_boundaries(self):
        for border in self.mesh_parameters.borders:
            if border.boundary_type == BoundaryType.Dirichlet:
                self.process_boundary_condition(border, self.dirichlet)
            elif border.boundary_type == BoundaryType.Neumann:
                self.process_boundary_condition(border, self.neumann)
            elif border.boundary_type == BoundaryType.Newton:
                self.process_boundary_condition(border, self.newton)
            elif border.boundary_type == BoundaryType.Null:
                continue
            else:
                self.process_boundary_condition(border, self.dirichlet)

    def process_boundary_condition(self, border: Border, conditions: list):
        x_splits = self.mesh_parameters.abscissa_splits
        y_splits = self.mesh_parameters.ordinate_splits

        nx = self.mesh_parameters.abscissa_points_count
        total_nx = self.mesh_parameters.abscissa_splits + 1

        ys = border.points_indices[0] // nx
        xs = border.points_indices[0] - ys * nx
        ye = border.points_indices[1] // nx
        xe = border.points_indices[1] - ye * nx

        xs = self.ir[xs]
        xe = self.ir[xe]
        ys = self.iz[ys]
        ye = self.iz[ye]

        (xs, xe) = (xe, xs) if xe < xs else (xs, xe)
        (ys, ye) = (ye, ys) if ye < ys else (ys, ye)

        # Horizontal line (начало и конец по y равны)
        if ys == ye:
            # ys=ye=0 значит это нижняя граница (т.к. индекс первой точки по оси Y равен 0)
            local_border_index = 0 if ys == 0 else 3
            ielem = 0 if local_border_index == 0 else x_splits * (y_splits - 1)

            for j in range(len(self.points)):
                iy = j // total_nx
                ix = j - iy * total_nx

                # if point is not on the line
                if iy != ys or ix < xs or ix > xe:
                    continue

                formula = self.mesh_parameters.boundary_formulas[border.formula_index]

                if isinstance(formula, BoundaryFormulaS3):
                    conditions.append(BoundaryConditionS3(ielem, local_border_index, formula.value, formula.beta))
                    ielem += 1
                else:
                    conditions.append(BoundaryCondition(ielem, local_border_index, formula.value))
                    ielem += 1

                # идем по всем точкам горизонтальной оси
                while ix + 1 != xe:
                    j += 1
                    ix = j - iy * total_nx

                    formula = self.mesh_parameters.boundary_formulas[border.formula_index]

                    if isinstance(formula, BoundaryFormulaS3):
                        conditions.append(BoundaryConditionS3(ielem, local_border_index, formula.value, formula.beta))
                        ielem += 1
                    else:
                        conditions.append(BoundaryCondition(ielem, local_border_index, formula.value))
                        ielem += 1

                break
        # Vertical line
        else:
            # xs=xe=0 значит это левая граница (т.к. индекс первой точки по оси X равен 0)
            local_border_index = 1 if xs == 0 else 2
            ielem = 0 if local_border_index == 1 else x_splits - 1

            for j in range(len(self.points)):
                iy = j // total_nx
                ix = j - iy * total_nx

                # If point is not on the line
                if ix != xs or iy < ys or iy > ye:
                    continue

                formula = self.mesh_parameters.boundary_formulas[border.formula_index]

                if isinstance(formula, BoundaryFormulaS3):
                    conditions.append(BoundaryConditionS3(ielem, local_border_index, formula.value, formula.beta))
                else:
                    conditions.append(BoundaryCondition(ielem, local_border_index, formula.value))

                # идем по всем точкам вертикальной оси
                while iy + 1 != ye:
                    ielem += x_splits
                    j += total_nx
                    iy = j // total_nx

                    formula = self.mesh_parameters.boundary_formulas[border.formula_index]

                    if isinstance(formula, BoundaryFormulaS3):
                        conditions.append(BoundaryConditionS3(ielem, local_border_index, formula.value, formula.beta))
                    else:
                        conditions.append(BoundaryCondition(ielem, local_border_index, formula.value))

                break

    def get_mesh(self):
        return Mesh(
            self.points,
            self.elements,
            self.mesh_parameters.area_properties,
            self.dirichlet,
            self.neumann,
            self.newton
        )