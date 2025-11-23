from mesh.mesh import Mesh

class Utils:
    @staticmethod
    def save_mesh(mesh: Mesh):
        # points
        with open("output/points", "w") as file:
            for p in mesh.points:
                file.write(f"{p.r} {p.z}\n")

        # elements
        with open("output/elements", "w") as file:
            for e in mesh.elements:
                nodes = e.physical_nodes_indices
                file.write(f"{nodes[0]} {nodes[1]} {nodes[2]} {nodes[3]}\n")

        # dirichlet
        with open("output/dirichlet", "w") as file:
            for d in mesh.dirichlet:
                file.write(f"{d.element} {d.local_border} {d.value}\n")

        # neumann
        with open("output/neumann", "w") as file:
            for n in mesh.neumann:
                file.write(f"{n.element} {n.local_border} {n.value}\n")

        #newton
        with open("output/newton", "w") as file:
            for n in mesh.newton:
                file.write(f"{n.element} {n.local_border} {n.value}\n")

    @staticmethod
    def save_basis_info(mesh: Mesh):
        with open("output/basis", "w") as file:
            for e in mesh.elements:
                nodes = e.basis_indices
                file.write(" ".join(str(node) for node in nodes) + "\n")

    @staticmethod
    def print_vector(vector: list, path: str):
        with open(f"output/{path}", "w") as file:
            for i in range(len(vector)):
                file.write(f"{vector[i]}\n")