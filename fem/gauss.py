import math

class Gauss:
    points = [
        0.0,
        1.0 / 3.0 * math.sqrt(5 - 2 * math.sqrt(10.0 / 7.0)),
        -1.0 / 3.0 * math.sqrt(5 - 2 * math.sqrt(10.0 / 7.0)),
        1.0 / 3.0 * math.sqrt(5 + 2 * math.sqrt(10.0 / 7.0)),
        -1.0 / 3.0 * math.sqrt(5 + 2 * math.sqrt(10.0 / 7.0))
    ]

    weights = [
        128.0 / 225.0,
        (322.0 + 13.0 * math.sqrt(70.0)) / 900.0,
        (322.0 + 13.0 * math.sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * math.sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * math.sqrt(70.0)) / 900.0
    ]