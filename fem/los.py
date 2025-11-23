from fem.sparse_matrix import SparseMatrix
from mesh.point import Point

class Los:
    def __init__(self, max_iterations: int, eps: float):
        self.max_iterations = max_iterations
        self.eps = eps
        self.solution: list[float] = []
        self.iterations_count: int = 0

    def compute(self, matrix: SparseMatrix, right_part: list[float]):
        try:
            n = len(right_part)
            self.solution = [0.0] * n

            z = [0.0] * n
            r = [0.0] * n
            p = [0.0] * n
            product = [0.0] * n

            # r = right_part - A * solution
            SparseMatrix.dot(matrix, self.solution, product)
            for i in range(n):
                r[i] = right_part[i] - product[i]

            z = r.copy()
            SparseMatrix.dot(matrix, z, p)

            square_norm = Point.dot(r, r)

            for self.iterations_count in range(self.max_iterations):
                if square_norm < self.eps:
                    break

                alpha = Point.dot(p, r) / Point.dot(p, p)

                for i in range(n):
                    self.solution[i] += alpha * z[i]
                    r[i] -= alpha * p[i]

                square_norm = Point.dot(r, r)

                if square_norm < self.eps:
                    break

                SparseMatrix.dot(matrix, r, product)
                beta = -Point.dot(p, product) / Point.dot(p, p)

                for i in range(n):
                    z[i] = r[i] + beta * z[i]
                    p[i] = product[i] + beta * p[i]

        except Exception as e:
            print(f"We had problem: {e}")
            raise