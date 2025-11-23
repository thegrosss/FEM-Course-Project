class SparseMatrix:
    def __init__(self, ig: list[int], jg: list[int]):
        self.ig = ig
        self.jg = jg
        self.di = [0.0] * (len(ig) - 1)
        self.gg = [0.0] * len(jg)
        self.size = len(self.di)

    def add(self, i: int, j: int, value: float):
        if i == j:
            self.di[i] += value
        elif i > j:
            for idx in range(self.ig[i], self.ig[i + 1]):
                if self.jg[idx] != j: continue
                self.gg[idx] += value

    def dot(self, vector: list[float], product: list[float] = None):
        if self.size != len(vector):
            raise Exception("Size of matrix not equal to size of vector")

        if product is None:
            product = [0.0] * len(vector)
        else:
            for i in range(len(product)):
                product[i] = 0.0

        for i in range(self.size):
            product[i] += self.di[i] * vector[i]

            for j in range(self.ig[i], self.ig[i + 1]):
                j_index = self.jg[j]
                product[i] += self.gg[j] * vector[j_index]
                product[j_index] += self.gg[j] * vector[i]

        return product

    def print_dense(self, path: str):
        a = [[0.0 for _ in range(self.size)] for _ in range(self.size)]

        for i in range(self.size):
            a[i][i] = self.di[i]

            for idx in range(self.ig[i], self.ig[i + 1]):
                j = self.jg[idx]
                a[i][j] = self.gg[idx]
                a[j][i] = self.gg[idx]

        with open(path, "w") as file:
            for i in range(self.size):
                for j in range(self.size):
                    file.write(f"{a[i][j]:.7f}\t")
                file.write("\n")

    def clear(self):
        for i in range(self.size):
            self.di[i] = 0.0
        for i in range(len(self.gg)):
            self.gg[i] = 0.0