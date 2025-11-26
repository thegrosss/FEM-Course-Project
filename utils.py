from mesh.mesh import Mesh
from mesh.point import Point

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
    def save_solution(mesh: Mesh, solution: list[float]):
        dict_solution: dict[int, tuple[Point, float]] = {}

        for element in mesh.elements:
            for i in range(9):
                g = element.get_global_basis_index(i)

                if g in dict_solution:
                    continue

                global_point = element.get_basis_node_position(i, lambda idx: mesh.points[idx])
                dict_solution[g] = (global_point, solution[g])

        ordered_solution = sorted(dict_solution.items(), key=lambda x: x[0])

        with open("output/solution", "w") as file:
            for global_idx, (point, value) in ordered_solution:
                file.write(f"{point.r} {point.z} {value}\n")

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