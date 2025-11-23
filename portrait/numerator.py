from mesh.mesh import Mesh

class Numerator:

    @staticmethod
    def numerate_basis_functions(mesh: Mesh):
        nx = mesh.elements[0].physical_nodes_indices[2]

        for ielem in range(len(mesh.elements)):
            element = mesh.elements[ielem]
            k = 2 * (ielem // (nx - 1)) * (2 * nx - 1) + 2 * (ielem % (nx - 1))

            element.set_basis_index(0, k)
            element.set_basis_index(1, k + 1)
            element.set_basis_index(2, k + 2)
            element.set_basis_index(3, k + 2 * nx - 1)
            element.set_basis_index(4, k + 2 * nx)
            element.set_basis_index(5, k + 2 * nx + 1)
            element.set_basis_index(6, k + 4 * nx - 2)
            element.set_basis_index(7, k + 4 * nx - 1)
            element.set_basis_index(8, k + 4 * nx)
