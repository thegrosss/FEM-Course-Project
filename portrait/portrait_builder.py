from mesh.mesh import Mesh

class PortraitBuilder:
    @staticmethod
    def generate_portrait(mesh: Mesh) -> tuple[list[int], list[int]]:
        func_count = mesh.elements[-1].basis_indices[-1] + 1
        connectivity_list: list[list[int]] = []

        for i in range(func_count):
            connectivity_list.append([])

        for element in mesh.elements:
            for node_to_insert in element.basis_indices:
                for pos_to_insert in element.basis_indices:
                    if node_to_insert < pos_to_insert:
                        if node_to_insert not in connectivity_list[pos_to_insert]:
                            connectivity_list[pos_to_insert].append(node_to_insert)

        for i in range(len(connectivity_list)):
            connectivity_list[i].sort()

        ig: list[int] = [0] * (len(connectivity_list) + 1)

        ig[0] = 0
        ig[1] = 0

        for i in range(1, len(connectivity_list)):
            ig[i + 1] = ig[i] + len(connectivity_list[i])

        jg: list[int] = [0] * ig[-1]

        j = 0
        for i in range(1, len(connectivity_list)):
            for it in connectivity_list[i]:
                jg[j] = it
                j += 1

        return ig, jg