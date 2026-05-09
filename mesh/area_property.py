class AreaProperty:
    def __init__(
        self,
        lmbda: float,
        gamma: float = 0.0,
        f = None,
        sigma: float = 0.0,
        hi: float = 0.0,
    ):
        self.lmbda = lmbda
        self.gamma = gamma
        self.sigma = sigma
        self.hi = hi
        self.f = f
