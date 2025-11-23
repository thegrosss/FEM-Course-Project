from mesh.edge import Edge
from mesh.point import Point

class BiquadraticQuadElement:
    basis_on_borders = [
        [0, 1, 2],  # низ
        [0, 3, 6],  # лево
        [2, 5, 8],  # право
        [6, 7, 8]  # верх
    ]

    def __init__(self, nodes: list[int], area_number: int):
        self.local_basis_to_physical_nodes: list[int] = [0, -1, 1, -1, -1, -1, 2, -1, 3]

        self.mid_edge_corner_pairs: list[tuple[int, int]] = [
            None,
            (0, 1),
            None,
            (0, 2),
            None,
            (1, 3),
            None,
            (2, 3),
            None
        ]

        self.area_number = area_number
        self.physical_nodes_indices = nodes
        self.basis_indices: list[int] = [0] * 9

        self.edges: list[Edge] = [
            Edge(nodes[0], nodes[1]),
            Edge(nodes[0], nodes[2]),
            Edge(nodes[1], nodes[3]),
            Edge(nodes[2], nodes[3])
        ]

    def get_global_node_index_for_basis(self, local_basis_index: int) -> int:
        if local_basis_index < 0 or local_basis_index > 8:
            raise IndexError("local_basis_index out of range")

        corner = self.local_basis_to_physical_nodes[local_basis_index]

        if corner == -1:
            return -1

        return self.physical_nodes_indices[corner]

    def set_basis_index(self, local_basis_index: int, global_index: int):
        self.basis_indices[local_basis_index] = global_index

    def get_global_basis_index(self, local_basis_index: int) -> int:
        return self.basis_indices[local_basis_index]

    def get_basis_node_position(self, local_basis_index: int, node_position_provider):
        corner = self.local_basis_to_physical_nodes[local_basis_index]

        if corner != -1:
            global_index = self.physical_nodes_indices[corner]
            return node_position_provider(global_index)

        pair = self.mid_edge_corner_pairs[local_basis_index]

        if pair is not None:
            a: Point = node_position_provider(self.physical_nodes_indices[pair[0]])
            b: Point = node_position_provider(self.physical_nodes_indices[pair[1]])

            return Point((a.r + b.r) * 0.5, (a.z + b.z) * 0.5)

        if local_basis_index == 4:
            p0: Point = node_position_provider(self.physical_nodes_indices[0])
            p1: Point = node_position_provider(self.physical_nodes_indices[1])
            p2: Point = node_position_provider(self.physical_nodes_indices[2])
            p3: Point = node_position_provider(self.physical_nodes_indices[3])

            return Point((p0.r + p1.r + p2.r + p3.r) * 0.25, (p0.z + p1.z + p2.z + p3.z) * 0.25)

        raise ValueError("Unhandled local basis index.")

    @classmethod
    def get_basis_by_border(cls, local_border_index: int):
        return cls.basis_on_borders[local_border_index]

    def __repr__(self):
        return f"BiquadraticQuadElement(nodes={self.physical_nodes_indices}, area={self.area_number})"