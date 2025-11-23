from fem.gauss import Gauss

class Integrator:
    @staticmethod
    def integration2D(f, rk: float, rk1: float, zk: float, zk1: float):
        hr = hz = 0
        s = 0

        for i in range(5):
            qi = Gauss.weights[i]
            hr = abs(rk1 - rk)
            pi = (rk + rk1 + Gauss.points[i] * hr) / 2.0

            for j in range(5):
                qj = Gauss.weights[j]
                hz = abs(zk1 - zk)
                pj = (zk + zk1 + Gauss.points[j] * hz) / 2.0

                s += qi * qj * f(pi, pj)

        return s * hr * hz / 4.0

    @staticmethod
    def integration1D(f, left_border: float, right_border: float):
        h = abs(right_border - left_border)
        s = 0

        for i in range(5):
            qi = Gauss.weights[i]
            pi = (left_border + right_border + Gauss.points[i] * h) / 2.0

            s += qi * f(pi)

        return s * h / 2.0