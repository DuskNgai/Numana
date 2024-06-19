import math
from pathlib import Path
import sys

import pytest
import sympy as sp

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numana.chapter5 as ch5

class TestNumericalDifferentiation(object):
    def output_numerical_differentiation(self, f: sp.Function, x: float, h: float):
        nd = ch5.NumericalDifferentiation(f)
        result = nd(x, h)
        print("The function is \033\13331mf(x) = {}\033\1330m and the point is \033\13331mx = {}\033\1330m with increment \033\13331mh = {}\033\1330m.".format(nd.symbol_f, x, h))
        print("Using two-points forward-difference formula for first-order derivative,     the result is: \033\13334m[{}]\033\1330m.".format(result[0]))
        print("Using two-points backward-difference formula for first-order derivative,    the result is: \033\13334m[{}]\033\1330m.".format(result[1]))
        print("Using three-points centered-difference formula for first-order derivative,  the result is: \033\13334m[{}]\033\1330m.".format(result[2]))
        print("Using five-points centered-difference formula for first-order derivative,   the result is: \033\13334m[{}]\033\1330m.".format(result[3]))
        print("Using symbolic differentiation formula for first-order derivative,          the result is: \033\13334m[{}]\033\1330m.".format(result[4]))
        print("Using three-points centered-difference formula for second-order derivative, the result is: \033\13334m[{}]\033\1330m.".format(result[5]))
        print("Using symbolic differentiation formula for second-order derivative,         the result is: \033\13334m[{}]\033\1330m.".format(result[6]))

    def test_numerical_differentiation(self):
        x = sp.symbols('x')
        self.output_numerical_differentiation(1.0 / x, 2.0, 0.1)
        # (1)
        self.output_numerical_differentiation(sp.log(x), 1.0, 1e-1)
        self.output_numerical_differentiation(sp.log(x), 1.0, 1e-2)
        self.output_numerical_differentiation(sp.log(x), 1.0, 1e-3)
        # (2)
        self.output_numerical_differentiation(sp.exp(x), 0.0, 1e-1)
        self.output_numerical_differentiation(sp.exp(x), 0.0, 1e-2)
        self.output_numerical_differentiation(sp.exp(x), 0.0, 1e-3)
        self.output_numerical_differentiation(sp.exp(x), 0.0, 1e-4)
        self.output_numerical_differentiation(sp.exp(x), 0.0, 1e-5)
        self.output_numerical_differentiation(sp.exp(x), 0.0, 1e-6)
        # (3), (4)
        self.output_numerical_differentiation(sp.sin(x), math.pi / 3, 1e-1)
        self.output_numerical_differentiation(sp.sin(x), math.pi / 3, 1e-2)
        self.output_numerical_differentiation(sp.sin(x), math.pi / 3, 1e-3)

class TestNewtonCotes(object):
    def output_newton_cotes(self, f: sp.Function, a: float, b: float):
        nc = ch5.NewtonCotes(f)
        result = nc(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(nc.symbol_f, a, b))
        for i in range(7):
            print("Using \033\13331m{0}\033\1330m order Newton-Cotes formula, the result is: \033\13334m[{1}]\033\1330m.".format(i + 1, result[i]))
        print("Using symbolic integration formula, the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))

    def test_newton_cotes(self):
        x = sp.symbols('x')
        self.output_newton_cotes(sp.log(x), 1.0, 2.0)
        self.output_newton_cotes(x ** 2,    0.0, 1.0)
        self.output_newton_cotes(sp.cos(x), 0.0, math.pi / 2.0)
        self.output_newton_cotes(sp.exp(x), 0.0, 1.0)

class TestCompositeNewtonCotes(object):
    def output_composite_newton_cotes(self, f: sp.Function, a: float, b: float, m: int):
        nc = ch5.CompositeNewtonCotes(f)
        result = nc(a, b, m)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m, we divide the it into \033\13331m[{}]\033\1330m intervals.".format(nc.symbol_f, a, b, m))
        print("Using \033\13331mComposite Trapezoid Rule\033\1330m,       the result is: \033\13334m[{}]\033\1330m.".format(result[0]))
        print("Using \033\13331mComposite Midpoint  Rule\033\1330m,       the result is: \033\13334m[{}]\033\1330m.".format(result[1]))
        print("Using \033\13331mComposite Three Midpoints Rule\033\1330m, the result is: \033\13334m[{}]\033\1330m.".format(result[2]))
        print("Using \033\13331mComposite Simpson's Rule\033\1330m,       the result is: \033\13334m[{}]\033\1330m.".format(result[3]))
        print("Using symbolic integration,                                the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))

    def test_composite_newton_cotes(self):
        x = sp.symbols('x')
        self.output_composite_newton_cotes(sp.log(x), 1.0, 2.0, 4)
        self.output_composite_newton_cotes(sp.sin(x) / x, 0.0, 1.0, 4)
        # (1), (2), (3)
        self.output_composite_newton_cotes(x ** 2, 0.0, 1.0, 1)
        self.output_composite_newton_cotes(x ** 2, 0.0, 1.0, 2)
        self.output_composite_newton_cotes(x ** 2, 0.0, 1.0, 4)
        self.output_composite_newton_cotes(sp.cos(x), 0.0, math.pi / 2.0, 1)
        self.output_composite_newton_cotes(sp.cos(x), 0.0, math.pi / 2.0, 2)
        self.output_composite_newton_cotes(sp.cos(x), 0.0, math.pi / 2.0, 4)
        self.output_composite_newton_cotes(sp.exp(x), 0.0, 1.0, 1)
        self.output_composite_newton_cotes(sp.exp(x), 0.0, 1.0, 2)
        self.output_composite_newton_cotes(sp.exp(x), 0.0, 1.0, 4)
        # (4)
        self.output_composite_newton_cotes(x * sp.exp(x), 0.0, 1.0, 1)
        self.output_composite_newton_cotes(x * sp.exp(x), 0.0, 1.0, 2)
        self.output_composite_newton_cotes(x * sp.exp(x), 0.0, 1.0, 4)
        # (5)
        self.output_composite_newton_cotes(1 / (1 + x ** 2), 0.0, 1.0, 1)
        self.output_composite_newton_cotes(1 / (1 + x ** 2), 0.0, 1.0, 2)
        self.output_composite_newton_cotes(1 / (1 + x ** 2), 0.0, 1.0, 4)
        # (6)
        self.output_composite_newton_cotes(x * sp.cos(x), 0.0, math.pi, 1)
        self.output_composite_newton_cotes(x * sp.cos(x), 0.0, math.pi, 2)
        self.output_composite_newton_cotes(x * sp.cos(x), 0.0, math.pi, 4)

        # (7)
        self.output_composite_newton_cotes(sp.sin(x) / x, 0.0, math.pi / 2.0, 16)
        self.output_composite_newton_cotes(sp.sin(x) / x, 0.0, math.pi / 2.0, 32)
        self.output_composite_newton_cotes((sp.exp(x) - 1) / sp.sin(x), 0.0, math.pi / 2.0, 16)
        self.output_composite_newton_cotes((sp.exp(x) - 1) / sp.sin(x), 0.0, math.pi / 2.0, 32)
        self.output_composite_newton_cotes(sp.atan(x) / x, 0.0, 1 / 2.0, 16)
        self.output_composite_newton_cotes(sp.atan(x) / x, 0.0, 1 / 2.0, 32)

class TestRomberg(object):
    def output_romberg(self, f: sp.Function, a: float, b: float, m: int):
        romberg = ch5.Romberg(f)
        result = romberg(a, b, m)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m, line of Romberg table is \033\13331m[{}]\033\1330m.".format(romberg.symbol_f, a, b, m))
        print("Using \033\13331mRomberg integration\033\1330m , the result is:")
        for i in range(m):
            for j in range(i + 1):
                print("\033\13334m[{:.16f}]\033\1330m".format(result[0][i][j]), end=" ")
            print()
        print("Using symbolic integration, the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))

    def test_romberg(self):
        x = sp.symbols('x')
        self.output_romberg(sp.log(x), 1.0, 2.0, 4)
        # (1)
        self.output_romberg(x ** 2, 0.0, 1.0, 3)
        self.output_romberg(sp.cos(x), 0.0, math.pi / 2.0, 3)
        self.output_romberg(sp.exp(x), 0.0, 1.0, 3)
        # (2)
        self.output_romberg(x * sp.exp(x), 0.0, 1.0, 3)
        self.output_romberg(1.0 / (1.0 + x ** 2), 0.0, 1.0, 3)
        self.output_romberg(x * sp.cos(x), 0.0, math.pi, 3)
        # (1)
        self.output_romberg(x / sp.sqrt(x ** 2 + 9), 0.0, 4.0, 5)
        self.output_romberg(x ** 3 / sp.sqrt(x ** 2 + 1), 0.0, 1.0, 5)
        self.output_romberg(x * sp.exp(x), 0.0, 1.0, 5)
        self.output_romberg(x ** 2 * sp.log(x), 1.0, 3.0, 5)
        self.output_romberg(x ** 2 * sp.sin(x), 0.0, math.pi, 5)
        self.output_romberg(x ** 3 / sp.sqrt(x ** 4 - 1), 2.0, 3.0, 5)
        self.output_romberg(1 / sp.sqrt(x ** 2 + 4), 0.0, 2 * math.sqrt(3), 5)
        self.output_romberg(x / sp.sqrt(x ** 4 + 1), 0.0, 1.0, 5)

class TestGaussLegendre(object):
    def output_gauss_legendre(self, f: sp.Function, a: float, b: float):
        gl = ch5.GaussLegendre(f)
        result = gl(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(gl.symbol_f, a, b))
        for i in range(5):
            print("Using \033\13331m{0}\033\1330m order Gauss-Legendre formula, the result is: \033\13334m[{1}]\033\1330m.".format(i + 1, result[i]))
        print("Using symbolic integration, the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))

    def test_gauss_legendre(self):
        x = sp.symbols('x')
        self.output_gauss_legendre(sp.exp(-x ** 2 / 2), -1.0, 1.0)
        self.output_gauss_legendre(sp.log(x), 1.0, 2.0)
        # (1) (2) (3)
        self.output_gauss_legendre((x ** 3 + 2 * x), -1.0, 1.0)
        self.output_gauss_legendre(x ** 4, -1.0, 1.0)
        self.output_gauss_legendre(sp.exp(x), -1.0, 1.0)
        self.output_gauss_legendre(sp.cos(sp.pi * x), -1.0, 1.0)
        # (4) (5)
        self.output_gauss_legendre(x / sp.sqrt(x ** 2 + 9), 0.0, 4.0)
        self.output_gauss_legendre(x ** 3 / sp.sqrt(x ** 2 + 1), 0.0, 1.0)
        self.output_gauss_legendre(x * sp.exp(x), 0.0, 1.0)
        self.output_gauss_legendre(x ** 2 * sp.log(x), 1.0, 3.0)
        # (6)
        self.output_gauss_legendre((x ** 3 + 2 * x), 0.0, 1.0)
        self.output_gauss_legendre(sp.log(x), 1.0, 4.0)
        self.output_gauss_legendre(x ** 5, -1.0, 2.0)
        self.output_gauss_legendre(sp.exp(-x ** 2 / 2), -3.0, 3.0)

if __name__ == "__main__":
    pytest.main(["-s", "test_ch5.py::TestNumericalDifferentiation::test_numerical_differentiation"])
    pytest.main(["-s", "test_ch5.py::TestNewtonCotes::test_newton_cotes"])
    pytest.main(["-s", "test_ch5.py::TestCompositeNewtonCotes::test_composite_newton_cotes"])
    pytest.main(["-s", "test_ch5.py::TestRomberg::test_romberg"])
    pytest.main(["-s", "test_ch5.py::TestGaussLegendre::test_gauss_legendre"])
