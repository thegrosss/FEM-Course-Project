class Edge:
    def __init__(self, node1: int, node2: int):
        self.node1 = node1
        self.node2 = node2

    def equal(self, other: "Edge"):
        return (self.node1 == other.node1 and self.node2 == other.node2
                or self.node1 == other.node2 and self.node2 == other.node1)