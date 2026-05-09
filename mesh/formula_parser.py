import math
from typing import Callable


SAFE_GLOBALS = {
    "__builtins__": {},
    "math": math,
    "exp": math.exp,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "log": math.log,
    "log10": math.log10,
    "sqrt": math.sqrt,
    "pi": math.pi,
    "e": math.e,
}


def parse_formula(formula: str) -> Callable[[float, float, float, float, float, float, float], float]:
    expression = formula.strip()

    if "=" in expression:
        expression = expression.split("=", 1)[1].strip()

    expression = expression.replace("^", "**")
    code = compile(expression, "<formula>", "eval")

    def value(
        x: float,
        y: float,
        t: float = 0.0,
        lmbda: float = 0.0,
        gamma: float = 0.0,
        sigma: float = 0.0,
        hi: float = 0.0,
    ) -> float:
        local_vars = {
            "x": x,
            "y": y,
            "r": x,
            "z": y,
            "t": t,
            "lmbda": lmbda,
            "lambda_": lmbda,
            "lam": lmbda,
            "gamma": gamma,
            "sigma": sigma,
            "hi": hi,
        }
        return float(eval(code, SAFE_GLOBALS, local_vars))

    return value
