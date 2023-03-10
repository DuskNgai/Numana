import math
from typing import Tuple

import numpy as np
import sympy as sp
import sympy.abc

class NumericalDifferentiation(object):
    """
    Numerical differentiation methods.

    Args:
        f (sp.Function): the function to differentiate.
    """

    def __init__(self, f: sp.Function):
        self.x = sympy.abc.x
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def __call__(self, x: float, h: float = 1e-6) -> Tuple[float, float, float, float, float, float]:
        """
        Evaluate the numerical differentiation methods at `x`.

        Args:
            x (float): The point to evaluate.
            h (float, optional): The step size. Defaults to 1e-6.

        Returns:
            Tuple[float, float, float, float, float, float]:
                The results of the two-point forward/backward method for first derivatives,
                the three-point centered method for first and second order derivatives,
                the five-point centered method for first derivatives,
                and the analytical solution.
        """

        return (
            self._twoPointForward(x, h),
            self._twoPointBackward(x, h),
            self._threePointCentered(x, h),
            self._fivePointCentered(x, h),
            sp.diff(self.symbol_f, self.x).subs(self.x, x),
            self._threePointCenteredSecondOrder(x, h),
            sp.diff(self.symbol_f, self.x, 2).subs(self.x, x),
        )

    def _twoPointForward(self, x: float, h: float) -> float:
        """Two-point forward-difference method."""
        return (self.numeric_f(x + h) - self.numeric_f(x)) / h

    def _twoPointBackward(self, x: float, h: float) -> float:
        """Two-point backward-difference method."""
        return (self.numeric_f(x) - self.numeric_f(x - h)) / h

    def _threePointCentered(self, x: float, h: float) -> float:
        """Three-point centered-difference method."""
        return (self.numeric_f(x + h) - self.numeric_f(x - h)) / (2 * h)

    def _fivePointCentered(self, x: float, h: float) -> float:
        """Five-point centered-difference method."""
        return (
            -self.numeric_f(x + 2 * h)
            + 8 * self.numeric_f(x + h)
            - 8 * self.numeric_f(x - h)
            + self.numeric_f(x - 2 * h)
        ) / (12 * h)

    def _threePointCenteredSecondOrder(self, x: float, h: float) -> float:
        """Three-point centered-difference method for second-order derivative."""
        return (self.numeric_f(x + h) - 2 * self.numeric_f(x) + self.numeric_f(x - h)) / (h * h)

class NewtonCotes(object):
    """
    Newton-Cotes methods for numerical integration.

    Args:
        f (sp.Function): the function to evaluate on some intervals.
    """

    def __init__(self, f: sp.Function):
        self.x = sympy.abc.x
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def __call__(self, a: float, b: float) -> Tuple[float, float, float, float, float, float, float, float]:
        """
        Evaluate the Newton-Cotes methods on the interval [a, b].

        Args:
            a (float): The lower bound of the interval.
            b (float): The upper bound of the interval.

        Returns:
            Tuple[float, float, float, float, float, float, float, float]: 
                The results of the first 7 order Newton-Cotes methods and the exact result.
        """

        assert not (math.isinf(a) or math.isnan(a)), "Invalid interval."
        assert not (math.isinf(b) or math.isnan(b)), "Invalid interval."

        return (
            self._firstOrder(a, b),
            self._secondOrder(a, b),
            self._thirdOrder(a, b),
            self._forthOrder(a, b),
            self._fifthOrder(a, b),
            self._sixthOrder(a, b),
            self._seventhOrder(a, b),
            sp.integrate(self.symbol_f, (self.x, a, b)),
        )

    def _firstOrder(self, a: float, b: float) -> float:
        """Trapezoid Rule."""
        h = (b - a)
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(b)
        return h * (y0 + y1) / 2.0

    def _secondOrder(self, a: float, b: float) -> float:
        """Simpson's Rule."""
        h = (b - a) / 2.0
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(a + h)
        y2 = self.numeric_f(b)
        return h * (y0 + 4.0 * y1 + y2) / 3.0

    def _thirdOrder(self, a: float, b: float) -> float:
        """Simpson's 3/8 Rule"""
        h = (b - a) / 3.0
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(a + h)
        y2 = self.numeric_f(b - h)
        y3 = self.numeric_f(b)
        return h * (y0 + 3.0 * (y1 + y2) + y3) * 3.0 / 8.0

    def _forthOrder(self, a: float, b: float) -> float:
        """Boole's Rule"""
        h = (b - a) / 4.0
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(a + h)
        y2 = self.numeric_f(a + 2.0 * h)
        y3 = self.numeric_f(b - h)
        y4 = self.numeric_f(b)
        return h * (7.0 * (y0 + y4) + 32.0 * (y1 + y3) + 12.0 * y2) * 2.0 / 45.0

    def _fifthOrder(self, a: float, b: float) -> float:
        h = (b - a) / 5.0
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(a + h)
        y2 = self.numeric_f(a + 2.0 * h)
        y3 = self.numeric_f(b - 2.0 * h)
        y4 = self.numeric_f(b - h)
        y5 = self.numeric_f(b)
        return h * (19.0 * (y0 + y5) + 75.0 * (y1 + y4) + 50.0 * (y2 + y3)) * 5.0 / 288.0

    def _sixthOrder(self, a: float, b: float) -> float:
        h = (b - a) / 6.0
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(a + h)
        y2 = self.numeric_f(a + 2.0 * h)
        y3 = self.numeric_f(a + 3.0 * h)
        y4 = self.numeric_f(b - 2.0 * h)
        y5 = self.numeric_f(b - h)
        y6 = self.numeric_f(b)
        return h * (41.0 * (y0 + y6) + 216.0 * (y1 + y5) + 27.0 * (y2 + y4) + 272.0 * y3) / 140.0

    def _seventhOrder(self, a: float, b: float) -> float:
        h = (b - a) / 7.0
        y0 = self.numeric_f(a)
        y1 = self.numeric_f(a + h)
        y2 = self.numeric_f(a + 2.0 * h)
        y3 = self.numeric_f(a + 3.0 * h)
        y4 = self.numeric_f(b - 3.0 * h)
        y5 = self.numeric_f(b - 2.0 * h)
        y6 = self.numeric_f(b - h)
        y7 = self.numeric_f(b)
        return h * (751.0 * (y0 + y7) + 3577.0 * (y1 + y6) + 1323.0 * (y2 + y5) + 2989.0 * (y3 + y4)) * 7.0 / 17280.0

class CompositeNewtonCotes(object):
    """
    Composite Newton-Cotes methods for numerical integration.

    Args:
        f (sp.Function): the function to evaluate on some intervals.
    """

    def __init__(self, f: sp.Function):
        self.x = sympy.abc.x
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def __call__(self, a: float, b: float, m: int) -> Tuple[float, float, float, float, float]:
        """
        Evaluate the composite Newton-Cotes methods on the interval [a, b], which is divided into m intervals.

        Args:
            a (float): The lower bound of the interval.
            b (float): The upper bound of the interval.
            m (int): The number of subintervals.

        Returns:
            Tuple[float, float, float, float, float]: 
                The results of the composite Newton-Cotes methods, in the order of trapezoid, midpoint, 3 midpoints, Simpson's, and the exact result.
        """

        assert not (math.isinf(a) or math.isnan(a)), "Invalid interval."
        assert not (math.isinf(b) or math.isnan(b)), "Invalid interval."
        assert m > 0, "Invalid number of intervals."

        return (
            self._trapezoid(a, b, m),
            self._midpoint(a, b, m),
            self._threeMidpoints(a, b, m),
            self._simpson(a, b, m),
            sp.integrate(self.symbol_f, (self.x, a, b)),
        )

    def _trapezoid(self, a: float, b: float, m: int) -> float:
        """Composite Trapezoid Rule."""
        h = (b - a) / m
        x = np.linspace(a, b, m + 1)
        y = self.numeric_f(x)
        return np.sum(y[:-1] + y[1:]) * h * 0.5

    def _midpoint(self, a: float, b: float, m: int) -> float:
        """Composite Midpoint Rule."""
        h = (b - a) / m
        h2 = h / 2.0
        x = np.linspace(a + h2, b - h2, m)
        y = self.numeric_f(x)
        return np.sum(y) * (b - a) / m

    def _threeMidpoints(self, a: float, b: float, m: int) -> float:
        """Composite Three Midpoints Rule."""
        h = (b - a) / m
        x = np.linspace(a, b, 4 * m + 1)
        y = self.numeric_f(x)
        return (2.0 * np.sum(y[1::2]) - np.sum(y[2::4])) * h / 3.0

    def _simpson(self, a: float, b: float, m: int) -> float:
        """Composite Simpson's Rule."""
        h = (b - a) / (2.0 * m)
        x = np.linspace(a, b, 2 * m + 1)
        y = self.numeric_f(x)
        # return (y[0] + y[-1] + 2.0 * np.sum(y[1::2]) + 2.0 * np.sum(y[1:-1])) * h / 3.0
        return (np.sum(y[:-2:2] + y[2::2]) + 4.0 * np.sum(y[1::2])) * h / 3.0

class Romberg(object):
    """
    Romberg methods for numerical integration.

    Args:
        f (sp.Function): the function to evaluate on some intervals.
    """

    def __init__(self, f: sp.Function):
        self.x = sympy.abc.x
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def __call__(self, a: float, b: float, m: int) -> Tuple[list[list[float]], float]:
        """
        Evaluate the m order Romberg methods on the interval [a, b].

        Args:
            a (float): The lower bound of the interval.
            b (float): The upper bound of the interval.
            m (int): The number of lines.

        Returns:
            Tuple[list[list[float]], float]: The results of all the m order Romberg method results and the exact result.
        """

        assert not (math.isinf(a) or math.isnan(a)), "Invalid interval."
        assert not (math.isinf(b) or math.isnan(b)), "Invalid interval."
        assert m > 0, "Invalid number of lines."

        h = (b - a) / 2.0
        y = self.numeric_f(np.linspace(a, b, 2 ** m + 1))
        R = [[] for _ in range(m)]

        R[0].append((y[0] + y[-1]) * h)
        for i in range(1, m):
            R[i].append((R[i - 1][0] / 2.0) + h * np.sum(y[2 ** (m - i)::2 ** (m - i + 1)]))
            for j in range(1, i + 1):
                R[i].append((4 ** j * R[i][j - 1] - R[i - 1][j - 1]) / (4 ** j - 1))
            h /= 2.0

        return (R, sp.integrate(self.symbol_f, (self.x, a, b)))

class GaussLegendre(object):
    """
    Gauss-Legendre methods for numerical integration.

    Args:
        f (sp.Function): the function to evaluate on some intervals.

    Attributes:
        ROOT (tuple[tuple[float]]): The roots of the first 5 Legendre polynomials.
        COEFFICIENT (tuple[tuple[float]]): The integration of the Lagrange interpolating polynomials of the roots on interval [-1, 1].
    """

    ROOT = (
        (0.0,),
        (-math.sqrt(3.0) / 3.0, math.sqrt(3.0) / 3.0),
        (-math.sqrt(15.0) / 5.0, 0.0, math.sqrt(15.0) / 5.0),
        (-math.sqrt(525.0 + 70.0 * math.sqrt(30.0)) / 35.0, -math.sqrt(525.0 - 70.0 * math.sqrt(30.0)) / 35.0, math.sqrt(525.0 - 70.0 * math.sqrt(30.0)) / 35.0, math.sqrt(525.0 + 70.0 * math.sqrt(30.0)) / 35.0),
        (-math.sqrt(245.0 + 14.0 * math.sqrt(70.0)) / 21.0, -math.sqrt(245.0 - 14.0 * math.sqrt(70.0)) / 21.0, 0.0, math.sqrt(245.0 - 14.0 * math.sqrt(70.0)) / 21.0, math.sqrt(245.0 + 14.0 * math.sqrt(70.0)) / 21.0)
    ) 

    COEFFICIENT = (
        (2.0,),
        (1.0, 1.0),
        (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0),
        (0.5 - math.sqrt(30.0) / 36.0, 0.5 + math.sqrt(30.0) / 36.0, 0.5 + math.sqrt(30.0) / 36.0, 0.5 - math.sqrt(30.0) / 36.0),
        ((322.0 - 13 * math.sqrt(70)) / 900.0, (322.0 + 13 * math.sqrt(70)) / 900.0, 128.0 / 225.0, (322.0 + 13 * math.sqrt(70)) / 900.0, (322.0 - 13 * math.sqrt(70)) / 900.0)
    )

    def __init__(self, f: sp.Function):
        self.x = sympy.abc.x
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def __call__(self, a: float, b: float) -> Tuple[float, float, float, float, float]:
        """
        Evaluate the Gauss-Legendre methods on the interval [a, b].

        Args:
            a (float): The lower bound of the interval.
            b (float): The upper bound of the interval.

        Returns:
            Tuple[float, float, float, float, float]: 
                The results of the first 4 order Gauss-Legendre methods and the exact result.
        """

        assert not (math.isinf(a) or math.isnan(a)), "Invalid interval."
        assert not (math.isinf(b) or math.isnan(b)), "Invalid interval."

        return (
            self._solver(a, b, 0),
            self._solver(a, b, 1),
            self._solver(a, b, 2),
            self._solver(a, b, 3),
            self._solver(a, b, 4),
            sp.integrate(self.symbol_f, (self.x, a, b)),
        )

    def _solver(self, a: float, b: float, order: int) -> float:
        result = 0.0
        for r, c in zip(self.ROOT[order], self.COEFFICIENT[order]):
            root = 0.5 * ((b - a) * r + b + a)
            result += c * self.numeric_f(root) * (b - a) * 0.5
        return result

class GaussChebyshev(object):
    """
    Gauss-Chebyshev methods for numerical integration.

    Args:
        f (sp.Function): the function to evaluate on some intervals.
    """

    def __init__(self, f: sp.Function):
        self.x = sympy.abc.x
        self.symbol_f = f
        self.numeric_f = sp.lambdify(self.x, f, "numpy")

    def __call__(self, a: float, b: float) -> Tuple[float, float, float, float, float]:
        """
        Evaluate the Gauss-Chebyshev methods on the interval [a, b].

        Args:
            a (float): The lower bound of the interval.
            b (float): The upper bound of the interval.

        Returns:
            Tuple[float, float, float, float, float]: 
                The results of the first 4 order Gauss-Chebyshev methods and the exact result.
        """

        assert not (math.isinf(a) or math.isnan(a)), "Invalid interval."
        assert not (math.isinf(b) or math.isnan(b)), "Invalid interval."

        return (
            self._solver(a, b, 0),
            self._solver(a, b, 1),
            self._solver(a, b, 2),
            self._solver(a, b, 3),
            self._solver(a, b, 4),
            sp.integrate(self.symbol_f, (self.x, a, b)),
        )

    def _solver(self, a: float, b: float, order: int) -> float:
        result = 0.0
        for i in range(1, order + 2):
            root = 0.5 * ((b - a) * math.cos((2 * i - 1) / (2 * (order + 1)) * math.pi) + b + a)
            result += self.numeric_f(root) * (b - a) * 0.5
        return result

