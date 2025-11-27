from fem.matrix import Matrix
from mesh.mesh import Mesh
from fem.basis import Basis
from fem.integrator import Integrator
from portrait.portrait_builder import PortraitBuilder
from fem.sparse_matrix import SparseMatrix
from mesh.biquadratic_quad_element import BiquadraticQuadElement
from typing import List, Tuple, Set

class MatrixAssembler:
    def __init__(self, mesh: Mesh):
        self.mesh = mesh

        self.ig, self.jg = PortraitBuilder.generate_portrait(mesh)

        self.G = Matrix(9, 9)    # stiffness (local)
        self.M = Matrix(9, 9)    # mass (local)
        self.local_b: List[float] = [0.0] * 9   # local RHS (M * f)
        self.local_f: List[float] = [0.0] * 9   # local source values

        self.global_b = [0.0] * (len(self.ig) - 1)
        self.global_matrix = SparseMatrix(self.ig, self.jg)

    def get_slae(self):
        self.assemble_global_slae()

        # сначала учитываются 2 и 3 краевые (порядок неважен)
        # 1-е краевые учитываются в последнюю очередь
        self.account_newton()
        self.account_neumann()
        self.account_dirichlet()

        return self.global_matrix, self.global_b

    def assemble_global_slae(self):
        self.global_matrix.clear()
        # обнуляем глобальную правую часть
        self.global_b = [0.0] * len(self.global_b)

        for ielem in range(len(self.mesh.elements)):
            mat = self.mesh.materials[self.mesh.elements[ielem].area_number]
            lmbda = mat.lmbda
            gamma = mat.gamma

            self.assemble_local_slae(ielem)

            for i in range(9):
                global_i = self.mesh.elements[ielem].get_global_basis_index(i)
                self.global_b[global_i] += self.local_b[i]

                for j in range(9):
                    global_j = self.mesh.elements[ielem].get_global_basis_index(j)
                    value = lmbda * self.G[i, j] + gamma * self.M[i, j]
                    self.global_matrix.add(global_i, global_j, value)

    def account_dirichlet(self):
        # сначала соберем все узлы для первого краевого в одном месте, чтобы проще учитывать
        all_dirichlet: List[Tuple[int, float]] = []
        processed_nodes: Set[int] = set()

        for d in self.mesh.dirichlet:
            element = self.mesh.elements[d.element]
            basis_by_border = BiquadraticQuadElement.get_basis_by_border(d.local_border)

            for local_basis_index in basis_by_border:
                global_basis = element.get_global_basis_index(local_basis_index)

                # Каждый узел нужно обработать только 1 раз
                if global_basis in processed_nodes:
                    continue
                processed_nodes.add(global_basis)

                global_point = element.get_basis_node_position(local_basis_index, lambda idx: self.mesh.points[idx])

                # на диагональ всегда ставим 1, а в правую часть ставим значение функции
                all_dirichlet.append((global_basis, d.value(global_point.r, global_point.z)))

        # если обратиться к последнему элементу и к его последней базисной функции, то можно узнать их количество
        # т.к. мы все пронумеровали последовательно
        f_count = self.mesh.elements[-1].basis_indices[-1] + 1
        bc1: List[int] = [-1 for _ in range(f_count)]

        for i in range(len(all_dirichlet)):
            bc1[all_dirichlet[i][0]] = i

        for i in range(f_count):
            if bc1[i] != -1:
                node, value = all_dirichlet[bc1[i]]

                self.global_matrix.di[node] = 1.0
                self.global_b[node] = value

                for j in range(self.global_matrix.ig[node], self.global_matrix.ig[node + 1]):
                    k = self.global_matrix.jg[j]

                    if bc1[k] == -1:
                        self.global_b[k] -= self.global_matrix.gg[j] * self.global_b[node]

                    self.global_matrix.gg[j] = 0.0
            else:
                for j in range(self.global_matrix.ig[i], self.global_matrix.ig[i + 1]):
                    k = self.global_matrix.jg[j]

                    if bc1[k] != -1:
                        self.global_b[i] -= self.global_matrix.gg[j] * self.global_b[k]
                        self.global_matrix.gg[j] = 0.0

    def account_neumann(self):
        # если 2х краевых нет, то и учитывать нечего
        if len(self.mesh.neumann) == 0:
            return

        for n in self.mesh.neumann:
            basis_by_border = BiquadraticQuadElement.get_basis_by_border(n.local_border)
            element = self.mesh.elements[n.element]
            border_start = element.get_basis_node_position(basis_by_border[0], lambda idx: self.mesh.points[idx])
            border_end = element.get_basis_node_position(basis_by_border[2], lambda idx: self.mesh.points[idx])
            rk = border_start.r
            rk1 = border_end.r
            zk = border_start.z
            zk1 = border_end.z

            # левая или правая граница
            if n.local_border == 1 or n.local_border == 2:
                for i in range(len(basis_by_border)):
                    local_basis = basis_by_border[i]
                    global_basis = element.get_global_basis_index(local_basis)
                    p = element.get_basis_node_position(local_basis, lambda idx: self.mesh.points[idx])
                    f = lambda z: Basis.psi_1d(i, zk, zk1, z)
                    self.global_b[global_basis] += rk * n.value(p.r, p.z) * Integrator.integration1D(f, zk, zk1)
            # нижняя или верхняя
            else:
                for i in range(len(basis_by_border)):
                    local_basis = basis_by_border[i]
                    global_basis = element.get_global_basis_index(local_basis)
                    p = element.get_basis_node_position(local_basis, lambda idx: self.mesh.points[idx])

                    # здесь еще в начале идет умножение на r, т.к это якобиан, при этом по горизонтальной оси изменяется r.
                    # поэтому нужно внести его под интеграл
                    f = lambda r: r * Basis.psi_1d(i, rk, rk1, r)
                    self.global_b[global_basis] += n.value(p.r, p.z) * Integrator.integration1D(f, rk, rk1)

    def account_newton(self):
        # если 3х краевых нет, то и учитывать нечего
        if len(self.mesh.newton) == 0:
            return

        total_local_mass_matrix = Matrix(3, 3)
        flow_vector = [0.0, 0.0, 0.0]
        local_vector = [0.0, 0.0, 0.0]

        for n in self.mesh.newton:
            basis_by_border = BiquadraticQuadElement.get_basis_by_border(n.local_border)
            element = self.mesh.elements[n.element]
            border_start = element.get_basis_node_position(basis_by_border[0], lambda idx: self.mesh.points[idx])
            border_end = element.get_basis_node_position(basis_by_border[2], lambda idx: self.mesh.points[idx])
            beta = n.beta
            rk = border_start.r
            rk1 = border_end.r
            zk = border_start.z
            zk1 = border_end.z

            for i in range(3):
                local_basis = basis_by_border[i]
                p = element.get_basis_node_position(local_basis, lambda idx: self.mesh.points[idx])
                flow_vector[i] = n.value(p.r, p.z)

            # левая или правая граница
            if n.local_border == 1 or n.local_border == 2:
                for i in range(3):
                    for j in range(i + 1):
                        f = lambda z: Basis.psi_1d(i, zk, zk1, z) * Basis.psi_1d(j, zk, zk1, z)
                        val = beta * rk * Integrator.integration1D(f, zk, zk1)
                        total_local_mass_matrix[i, j] = val
                        total_local_mass_matrix[j, i] = val
            # нижняя или верхняя
            else:
                for i in range(3):
                    for j in range(i + 1):
                        f = lambda r: r * Basis.psi_1d(i, rk, rk1, r) * Basis.psi_1d(j, rk, rk1, r)
                        val = beta * Integrator.integration1D(f, rk, rk1)
                        total_local_mass_matrix[i, j] = val
                        total_local_mass_matrix[j, i] = val

            local_vector = Matrix.dot(total_local_mass_matrix, flow_vector)

            for i in range(3):
                global_i = element.get_global_basis_index(basis_by_border[i])
                self.global_b[global_i] += local_vector[i]

                for j in range(i + 1):
                    global_j = element.get_global_basis_index(basis_by_border[j])
                    self.global_matrix.add(global_i, global_j, total_local_mass_matrix[i, j])

    def assemble_local_slae(self, ielem: int):
        element = self.mesh.elements[ielem]

        p0 = self.mesh.points[element.physical_nodes_indices[0]]
        p2 = self.mesh.points[element.physical_nodes_indices[-1]]

        rk = p0.r
        zk = p0.z
        rk1 = p2.r
        zk1 = p2.z

        # локальная матрица жесткости G
        for i in range(9):
            for j in range(i + 1):
                def f(r: float, z: float):
                    dpsi_i_r = Basis.d_psi(self.mesh, ielem, i, r, z, "r")
                    dpsi_j_r = Basis.d_psi(self.mesh, ielem, j, r, z, "r")
                    dpsi_i_z = Basis.d_psi(self.mesh, ielem, i, r, z, "z")
                    dpsi_j_z = Basis.d_psi(self.mesh, ielem, j, r, z, "z")
                    return (dpsi_i_r * dpsi_j_r + dpsi_i_z * dpsi_j_z) * r

                integral = Integrator.integration2D(f, rk, rk1, zk, zk1)
                self.G[i, j] = integral
                self.G[j, i] = integral

        # локальная матрица масс M
        for i in range(9):
            for j in range(i + 1):
                def f(r: float, z: float):
                    psi_i = Basis.psi(self.mesh, ielem, i, r, z)
                    psi_j = Basis.psi(self.mesh, ielem, j, r, z)
                    return psi_i * psi_j * r

                integral = Integrator.integration2D(f, rk, rk1, zk, zk1)
                self.M[i, j] = integral
                self.M[j, i] = integral

        f = self.mesh.materials[self.mesh.elements[ielem].area_number].f
        for i in range(9):
            point_i = element.get_basis_node_position(i, lambda idx: self.mesh.points[idx])
            self.local_f[i] = f(point_i.r, point_i.z)

        for i in range(9):
            self.local_b[i] = 0.0
            for j in range(9):
                self.local_b[i] += self.M[i, j] * self.local_f[j]
