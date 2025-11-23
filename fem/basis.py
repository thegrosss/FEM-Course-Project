from mesh.mesh import Mesh

class Basis:
    @staticmethod
    def psi_1d(function: int, rk: float, rk1: float, r: float):
        rkA = (rk + rk1) * 0.5
        hr2 = (rk1 - rk) * (rk1 - rk)

        if function == 0:
            return 2.0 / hr2 * (r - rkA) * (r - rk1)
        elif function == 1:
            return -4.0 / hr2 * (r - rk) * (r - rk1)
        elif function == 2:
            return 2.0 / hr2 * (r - rk) * (r - rkA)
        else:
            return 0.0

    @staticmethod
    def psi(mesh: Mesh, ielem: int, function: int, r: float, z: float):
        element = mesh.elements[ielem]

        p0 = mesh.points[element.physical_nodes_indices[0]]
        p2 = mesh.points[element.physical_nodes_indices[3]]

        rk = p0.r
        zk = p0.z
        rk1 = p2.r
        zk1 = p2.z

        rk_a = (rk + rk1) / 2.0
        zk_a = (zk + zk1) / 2.0

        hr2 = (rk1 - rk) * (rk1 - rk)
        hz2 = (zk1 - zk) * (zk1 - zk)

        if function == 0:
            return 2.0 / hr2 * (r - rk_a) * (r - rk1) * 2.0 / hz2 * (z - zk_a) * (z - zk1)
        elif function == 1:
            return (-4.0) / hr2 * (r - rk) * (r - rk1) * 2.0 / hz2 * (z - zk_a) * (z - zk1)
        elif function == 2:
            return 2.0 / hr2 * (r - rk) * (r - rk_a) * 2.0 / hz2 * (z - zk_a) * (z - zk1)
        elif function == 3:
            return 2.0 / hr2 * (r - rk_a) * (r - rk1) * (-4.0) / hz2 * (z - zk) * (z - zk1)
        elif function == 4:
            return (-4.0) / hr2 * (r - rk) * (r - rk1) * (-4.0) / hz2 * (z - zk) * (z - zk1)
        elif function == 5:
            return 2.0 / hr2 * (r - rk) * (r - rk_a) * (-4.0) / hz2 * (z - zk) * (z - zk1)
        elif function == 6:
            return 2.0 / hr2 * (r - rk_a) * (r - rk1) * 2.0 / hz2 * (z - zk) * (z - zk_a)
        elif function == 7:
            return (-4.0) / hr2 * (r - rk) * (r - rk1) * 2.0 / hz2 * (z - zk) * (z - zk_a)
        elif function == 8:
            return 2.0 / hr2 * (r - rk) * (r - rk_a) * 2.0 / hz2 * (z - zk) * (z - zk_a)
        else:
            return  0.0

    @staticmethod
    def d_psi(mesh: Mesh, ielem: int, function: int, r: float, z: float, var: str):
        element = mesh.elements[ielem]

        p0 = mesh.points[element.physical_nodes_indices[0]]
        p2 = mesh.points[element.physical_nodes_indices[3]]

        rk = p0.r
        zk = p0.z
        rk1 = p2.r
        zk1 = p2.z

        rk_a = (rk + rk1) / 2.0
        zk_a = (zk + zk1) / 2.0

        hr2 = (rk1 - rk) * (rk1 - rk)
        hz2 = (zk1 - zk) * (zk1 - zk)

        if var == "r":
            if function == 0:
                return 4 / (hr2 * hz2) * (z - zk_a) * (z - zk1) * ((r - rk_a) + (r - rk1))
            elif function == 1:
                return -8 / (hr2 * hz2) * (z - zk_a) * (z - zk1) * ((r - rk1) + (r - rk))
            elif function == 2:
                return 4 / (hr2 * hz2) * (z - zk_a) * (z - zk1) * ((r - rk_a) + (r - rk))
            elif function == 3:
                return -8 / (hr2 * hz2) * (z - zk) * (z - zk1) * ((r - rk_a) + (r - rk1))
            elif function == 4:
                return 16 / (hr2 * hz2) * (z - zk) * (z - zk1) * ((r - rk1) + (r - rk))
            elif function == 5:
                return -8 / (hr2 * hz2) * (z - zk) * (z - zk1) * ((r - rk_a) + (r - rk))
            elif function == 6:
                return 4 / (hr2 * hz2) * (z - zk) * (z - zk_a) * ((r - rk_a) + (r - rk1))
            elif function == 7:
                return -8 / (hr2 * hz2) * (z - zk) * (z - zk_a) * ((r - rk1) + (r - rk))
            elif function == 8:
                return 4 / (hr2 * hz2) * (z - zk) * (z - zk_a) * ((r - rk_a) + (r - rk))
            else:
                return 0.0
        elif var == "z":
            if function == 0:
                return 4 / (hr2 * hz2) * (r - rk_a) * (r - rk1) * ((z - zk_a) + (z - zk1))
            elif function == 1:
                return -8 / (hr2 * hz2) * (r - rk) * (r - rk1) * ((z - zk_a) + (z - zk1))
            elif function == 2:
                return 4 / (hr2 * hz2) * (r - rk) * (r - rk_a) * ((z - zk_a) + (z - zk1))
            elif function == 3:
                return -8 / (hr2 * hz2) * (r - rk_a) * (r - rk1) * ((z - zk1) + (z - zk))
            elif function == 4:
                return 16 / (hr2 * hz2) * (r - rk) * (r - rk1) * ((z - zk1) + (z - zk))
            elif function == 5:
                return -8 / (hr2 * hz2) * (r - rk) * (r - rk_a) * ((z - zk1) + (z - zk))
            elif function == 6:
                return 4 / (hr2 * hz2) * (r - rk_a) * (r - rk1) * ((z - zk_a) + (z - zk))
            elif function == 7:
                return -8 / (hr2 * hz2) * (r - rk) * (r - rk1) * ((z - zk_a) + (z - zk))
            elif function == 8:
                return 4 / (hr2 * hz2) * (r - rk) * (r - rk_a) * ((z - zk_a) + (z - zk))
            else:
                return 0.0
        else:
            return 0.0