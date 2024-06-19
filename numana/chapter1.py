import numpy as np
import sympy as sp

class Solver(object):
    epsilon: float = np.finfo(float).eps

    def __init__(self, f: sp.Function):
        self.x = sp.symbols("x")
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def symbolic(self) -> list:
        return sp.solve(self.symbol_f, self.x)

    def bisection(self, a: float, b: float) -> float:
        """
        Using bisection method to find a root in the interval `[a, b]`.
        Left endpoint is the anchor.

        Args:
            a (float): One evaluation point.
            b (float): One evaluation point.

        Returns:
            float: The root of the equation.
        """

        a, b = np.asarray(a), np.asarray(b)
        fa, fb = self.numeric_f(a), self.numeric_f(b)
        assert fa * fb < 0, "[{}, {}] is not a proper interval".format(a, b)
        tolerance = 2.0 * self.epsilon * max(1, abs(a))

        while abs(b - a) > tolerance:
            c = (a + b) / 2.0
            fc = self.numeric_f(c)
            if fc == 0.0:
                return c
            if fa * fc < 0:
                b = c
            else:
                a, fa = c, fc
        return (a + b) / 2

    def fixed_point_iteration(self, a: float) -> float:
        """
        Using fixed point iteration to find the root.

        Args:
            a (float): One evaluation point.

        Returns:
            float: The root of the equation.
        """

        a = np.asarray(a)
        b = self.numeric_f(a)
        tolerance = 2.0 * self.epsilon * max(1.0, abs(a))

        for _ in range(1000):
            if abs(b - a) > tolerance:
                a, b = b, self.numeric_f(b)
            else:
                break

        return b

    def newton(self, a: float) -> float:
        """
        Using Newton-Raffson's method to find the root.

        Args:
            a (float): One evaluation point.

        Returns:
            float: The root of the equation.
        """

        a, b = np.asarray(a), np.asarray(a)
        tolerance = 2.0 * self.epsilon * max(1.0, abs(a))
        df = sp.lambdify(self.x, sp.diff(self.symbol_f), "numpy")

        a, b = b, b - self.numeric_f(b) / df(b)
        while abs(b - a) > tolerance:
            a, b = b, b - self.numeric_f(b) / df(b)
        return a

    def secant(self, a: float, b: float) -> float:
        """
        Secant method (an improvement of Newton's method) for finding a root near `a` and `b`.

        Args:
            a (float): First evaluation point.
            b (float): Second evaluation point.

        Returns:
            float: The root of the equation.
        """

        a, b = np.asarray(a), np.asarray(b)
        fa, fb = self.numeric_f(a), self.numeric_f(b)
        tolerance = 2.0 * self.epsilon * max(1.0, abs(a))

        while abs(b - a) > tolerance:
            c = b - (fb * (b - a)) / (fb - fa)
            a, b = b, c
            fa, fb = fb, self.numeric_f(c)
        return b

    def regula_falsi(self, a: float, b: float) -> float:
        """
        A combination of bisection and secant method.

        Args:
            a (float): First evaluation point.
            b (float): Second evaluation point.

        Returns:
            float: The root of the equation.
        """

        a, b = np.asarray(a), np.asarray(b)
        fa, fb = self.numeric_f(a), self.numeric_f(b)
        tolerance = 2.0 * self.epsilon * max(1.0, abs(a))

        while abs(b - a) > tolerance:
            c = b - (fb * (b - a)) / (fb - fa)
            fc = self.numeric_f(c)
            if fc == 0.0:
                return c
            if fa * fc < 0.0:
                b, fb = c, fc
            else:
                a, fa = c, fc
        return b if abs(fb) < abs(self.numeric_f((a + b) / 2.0)) else (a + b) / 2.0

    def inverse_interpolation(self, a: float, b: float, c: float) -> float:
        """
        Inverse quadratic interpolation method for solving equations.

        Args:
            a (float): First evaluation point.
            b (float): Second evaluation point.
            c (float): Third evaluation point.

        Returns:
            float: The root of the equation.
        """

        a, b, c = np.asarray(a), np.asarray(b), np.asarray(c)
        fa, fb, fc = self.numeric_f(a), self.numeric_f(b), self.numeric_f(c)
        tolerance = 2.0 * self.epsilon * max(1.0, abs(a))

        while abs(b - c) > tolerance:
            AB = fa / fb
            AC = fa / fc
            BA = fb / fa
            BC = fb / fc
            CA = fc / fa
            CB = fc / fb

            d = -(a * (BA - CA) + b * (CB - AB) + c * (AC - BC)) / ((AB - 1) * (BC - 1) * (CA - 1))
            a, b, c = b, c, d
            fa, fb, fc = fb, fc, self.numeric_f(d)

        return c
