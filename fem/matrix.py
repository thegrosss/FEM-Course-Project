class Matrix:
    def __init__(self, n_rows: int, n_columns: int):
        self.rows = n_rows
        self.columns = n_columns
        self.storage = [[0.0] * n_columns for _ in range(n_rows)]

    def __getitem__(self, index):
        i, j = index
        return self.storage[i][j]

    def __setitem__(self, index, value):
        i, j = index
        self.storage[i][j] = value

    def fill(self, value: float):
        for i in range(self.rows):
            for j in range(self.columns):
                self.storage[i][j] = value

    def dot(self, vector: list[float], product: list[float]):
        if self.columns != len(vector):
            raise Exception("Numbers of columns not equal to size of vector")

        product = [0.0 for _ in range(len(vector))]
        for i in range(self.rows):
            for j in range(self.columns):
                product[i] += self.storage[i][j] * vector[j]