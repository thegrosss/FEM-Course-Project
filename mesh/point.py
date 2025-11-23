class Point:
    def __init__(self, r: float, z: float):
        self.r = r
        self.z = z

    @staticmethod
    def dot(a: list[float], b: list[float]):
        return sum(a_i * b_i for a_i, b_i in zip(a, b))