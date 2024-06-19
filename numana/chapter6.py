from typing import Tuple

import numpy as np
import sympy as sp

class EulerMethod(object):
    """
    Integrator for solving `y' = f(x, y)` with initial value `y(0)`.

    Args:
        f: The function `f(x, y)` in the equation.
    """
    def __init__(self, f: sp.Function):
        self.x = sp.Symbol('x')
        self.y = sp.Symbol('y')
        self.symbol_f = f
        self.numeric_f = sp.lambdify([self.x, self.y], f, "numpy")

    def explicit(self, initial: float, num_steps: int = 100, interval: Tuple[float, float] = (0.0, 1.0)) -> np.ndarray:
        """
        Using explicit Euler method.

        Args:
            initial (float): the initial value `y(0)`.
            num_steps (int): the number of steps.
            interval (Tuple[float, float]): the interval for solving.

        Returns:
            (np.ndarray): Solution at each step.
        """
        initial = np.asarray(initial)
        assert isinstance(num_steps, int) and num_steps > 0, "Number of steps must be a positive integer."
        start, end = interval
        assert isinstance(start, (int, float)) and isinstance(end, (int, float)) and start < end, "Not a valid interval."

        h = (end - start) / num_steps
        xs = np.linspace(start, end, num_steps + 1)
        ys = np.empty(num_steps + 1)
        ys[0] = initial
        for i in range(num_steps):
            ys[i + 1] = ys[i] + h * self.numeric_f(xs[i], ys[i])

        return ys

    def trapezoid(self, initial: float, num_steps: int = 100, interval: Tuple[float, float] = (0.0, 1.0)) -> np.ndarray:
        """
        Using trapezoid Euler method.
        """
        initial = np.asarray(initial)
        assert isinstance(num_steps, int) and num_steps > 0, "Number of steps must be a positive integer."
        start, end = interval
        assert isinstance(start, (int, float)) and isinstance(end, (int, float)) and start < end, "Not a valid interval."

        h = (end - start) / num_steps
        xs = np.linspace(start, end, num_steps + 1)
        ys = np.empty(num_steps + 1)
        ys[0] = initial
        for i in range(num_steps):
            mid = self.numeric_f(xs[i], ys[i])
            ys[i + 1] = ys[i] + 0.5 * h * (mid + self.numeric_f(xs[i] + h, ys[i] + h * mid))

        return ys
