import math
import os
import sys
sys.path.append(os.pardir)

import pytest
import sympy as sp
import sympy.abc

import numana.chapter5 as ch5

class TestNumericalDifferentiation(object):
    def outputNumericalDifferentiation(self, f: sp.Function, x: float, h: float):
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

    def testNumericalDifferentiation(self):
        x = sympy.abc.x
        self.outputNumericalDifferentiation(1.0 / x, 2.0, 0.1)
        # (1)
        self.outputNumericalDifferentiation(sp.log(x), 1.0, 1e-1)
        self.outputNumericalDifferentiation(sp.log(x), 1.0, 1e-2)
        self.outputNumericalDifferentiation(sp.log(x), 1.0, 1e-3)
        # (2)
        self.outputNumericalDifferentiation(sp.exp(x), 0.0, 1e-1)
        self.outputNumericalDifferentiation(sp.exp(x), 0.0, 1e-2)
        self.outputNumericalDifferentiation(sp.exp(x), 0.0, 1e-3)
        self.outputNumericalDifferentiation(sp.exp(x), 0.0, 1e-4)
        self.outputNumericalDifferentiation(sp.exp(x), 0.0, 1e-5)
        self.outputNumericalDifferentiation(sp.exp(x), 0.0, 1e-6)
        # (3), (4)
        self.outputNumericalDifferentiation(sp.sin(x), math.pi / 3, 1e-1)
        self.outputNumericalDifferentiation(sp.sin(x), math.pi / 3, 1e-2)
        self.outputNumericalDifferentiation(sp.sin(x), math.pi / 3, 1e-3)

class TestNewtonCotes(object):
    def outputNewtonCotes(self, f: sp.Function, a: float, b: float):
        nc = ch5.NewtonCotes(f)
        result = nc(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(nc.symbol_f, a, b))
        for i in range(7):
            print("Using \033\13331m{0}\033\1330m order Newton-Cotes formula, the result is: \033\13334m[{1}]\033\1330m.".format(i + 1, result[i]))
        print("Using symbolic integration formula, the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))
  
    def testNewtonCotes(self):
        x = sympy.abc.x
        self.outputNewtonCotes(sp.log(x), 1.0, 2.0)
        self.outputNewtonCotes(x ** 2,    0.0, 1.0)
        self.outputNewtonCotes(sp.cos(x), 0.0, math.pi / 2.0)
        self.outputNewtonCotes(sp.exp(x), 0.0, 1.0)

class TestCompositeNewtonCotes(object):
    def outputCompositeNewtonCotes(self, f: sp.Function, a: float, b: float, m: int):
        nc = ch5.CompositeNewtonCotes(f)
        result = nc(a, b, m)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m, we divide the it into \033\13331m[{}]\033\1330m intervals.".format(nc.symbol_f, a, b, m))
        print("Using \033\13331mComposite Trapezoid Rule\033\1330m,       the result is: \033\13334m[{}]\033\1330m.".format(result[0]))
        print("Using \033\13331mComposite Midpoint  Rule\033\1330m,       the result is: \033\13334m[{}]\033\1330m.".format(result[1]))
        print("Using \033\13331mComposite Three Midpoints Rule\033\1330m, the result is: \033\13334m[{}]\033\1330m.".format(result[2]))
        print("Using \033\13331mComposite Simpson's Rule\033\1330m,       the result is: \033\13334m[{}]\033\1330m.".format(result[3]))
        print("Using symbolic integration,                                the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))
  
    def testCompositeNewtonCotes(self):
        x = sympy.abc.x
        self.outputCompositeNewtonCotes(sp.log(x), 1.0, 2.0, 4)
        self.outputCompositeNewtonCotes(sp.sin(x) / x, 0.0, 1.0, 4)
        # (1), (2), (3)
        self.outputCompositeNewtonCotes(x ** 2, 0.0, 1.0, 1)
        self.outputCompositeNewtonCotes(x ** 2, 0.0, 1.0, 2)
        self.outputCompositeNewtonCotes(x ** 2, 0.0, 1.0, 4)
        self.outputCompositeNewtonCotes(sp.cos(x), 0.0, math.pi / 2.0, 1)
        self.outputCompositeNewtonCotes(sp.cos(x), 0.0, math.pi / 2.0, 2)
        self.outputCompositeNewtonCotes(sp.cos(x), 0.0, math.pi / 2.0, 4)
        self.outputCompositeNewtonCotes(sp.exp(x), 0.0, 1.0, 1)
        self.outputCompositeNewtonCotes(sp.exp(x), 0.0, 1.0, 2)
        self.outputCompositeNewtonCotes(sp.exp(x), 0.0, 1.0, 4)
        # (4)
        self.outputCompositeNewtonCotes(x * sp.exp(x), 0.0, 1.0, 1)
        self.outputCompositeNewtonCotes(x * sp.exp(x), 0.0, 1.0, 2)
        self.outputCompositeNewtonCotes(x * sp.exp(x), 0.0, 1.0, 4)
        # (5)
        self.outputCompositeNewtonCotes(1 / (1 + x ** 2), 0.0, 1.0, 1)
        self.outputCompositeNewtonCotes(1 / (1 + x ** 2), 0.0, 1.0, 2)
        self.outputCompositeNewtonCotes(1 / (1 + x ** 2), 0.0, 1.0, 4)
        # (6)
        self.outputCompositeNewtonCotes(x * sp.cos(x), 0.0, math.pi, 1)
        self.outputCompositeNewtonCotes(x * sp.cos(x), 0.0, math.pi, 2)
        self.outputCompositeNewtonCotes(x * sp.cos(x), 0.0, math.pi, 4)

        # (7)
        self.outputCompositeNewtonCotes(sp.sin(x) / x, 0.0, math.pi / 2.0, 16)
        self.outputCompositeNewtonCotes(sp.sin(x) / x, 0.0, math.pi / 2.0, 32)
        self.outputCompositeNewtonCotes((sp.exp(x) - 1) / sp.sin(x), 0.0, math.pi / 2.0, 16)
        self.outputCompositeNewtonCotes((sp.exp(x) - 1) / sp.sin(x), 0.0, math.pi / 2.0, 32)
        self.outputCompositeNewtonCotes(sp.atan(x) / x, 0.0, 1 / 2.0, 16)
        self.outputCompositeNewtonCotes(sp.atan(x) / x, 0.0, 1 / 2.0, 32)

class TestRomberg(object):
    def outputRomberg(self, f: sp.Function, a: float, b: float, m: int):
        romberg = ch5.Romberg(f)
        result = romberg(a, b, m)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m, line of Romberg table is \033\13331m[{}]\033\1330m.".format(romberg.symbol_f, a, b, m))
        print("Using \033\13331mRomberg integration\033\1330m , the result is:")
        for i in range(m):
            for j in range(i + 1):
                print("\033\13334m[{:.16f}]\033\1330m".format(result[0][i][j]), end=" ")
            print()
        print("Using symbolic integration, the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))

    def testRomberg(self):
        x = sympy.abc.x
        self.outputRomberg(sp.log(x), 1.0, 2.0, 4)
        # (1)
        self.outputRomberg(x ** 2, 0.0, 1.0, 3)
        self.outputRomberg(sp.cos(x), 0.0, math.pi / 2.0, 3)
        self.outputRomberg(sp.exp(x), 0.0, 1.0, 3)
        # (2)
        self.outputRomberg(x * sp.exp(x), 0.0, 1.0, 3)
        self.outputRomberg(1.0 / (1.0 + x ** 2), 0.0, 1.0, 3)
        self.outputRomberg(x * sp.cos(x), 0.0, math.pi, 3)
        # (1)
        self.outputRomberg(x / sp.sqrt(x ** 2 + 9), 0.0, 4.0, 5)
        self.outputRomberg(x ** 3 / sp.sqrt(x ** 2 + 1), 0.0, 1.0, 5)
        self.outputRomberg(x * sp.exp(x), 0.0, 1.0, 5)
        self.outputRomberg(x ** 2 * sp.log(x), 1.0, 3.0, 5)
        self.outputRomberg(x ** 2 * sp.sin(x), 0.0, math.pi, 5)
        self.outputRomberg(x ** 3 / sp.sqrt(x ** 4 - 1), 2.0, 3.0, 5)
        self.outputRomberg(1 / sp.sqrt(x ** 2 + 4), 0.0, 2 * math.sqrt(3), 5)
        self.outputRomberg(x / sp.sqrt(x ** 4 + 1), 0.0, 1.0, 5)

class TestGaussLegendre(object):
    def outputGaussLegendre(self, f: sp.Function, a: float, b: float):
        gl = ch5.GaussLegendre(f)
        result = gl(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(gl.symbol_f, a, b))
        for i in range(5):
            print("Using \033\13331m{0}\033\1330m order Gauss-Legendre formula, the result is: \033\13334m[{1}]\033\1330m.".format(i + 1, result[i]))
        print("Using symbolic integration, the result is: \033\13334m[{}]\033\1330m.".format(result[-1]))

    def testGaussLegendre(self):
        x = sympy.abc.x
        self.outputGaussLegendre(sp.exp(-x ** 2 / 2), -1.0, 1.0)
        self.outputGaussLegendre(sp.log(x), 1.0, 2.0)
        # (1) (2) (3)
        self.outputGaussLegendre((x ** 3 + 2 * x), -1.0, 1.0)
        self.outputGaussLegendre(x ** 4, -1.0, 1.0)
        self.outputGaussLegendre(sp.exp(x), -1.0, 1.0)
        self.outputGaussLegendre(sp.cos(sp.pi * x), -1.0, 1.0)
        # (4) (5)
        self.outputGaussLegendre(x / sp.sqrt(x ** 2 + 9), 0.0, 4.0)
        self.outputGaussLegendre(x ** 3 / sp.sqrt(x ** 2 + 1), 0.0, 1.0)
        self.outputGaussLegendre(x * sp.exp(x), 0.0, 1.0)
        self.outputGaussLegendre(x ** 2 * sp.log(x), 1.0, 3.0)
        # (6)
        self.outputGaussLegendre((x ** 3 + 2 * x), 0.0, 1.0)
        self.outputGaussLegendre(sp.log(x), 1.0, 4.0)
        self.outputGaussLegendre(x ** 5, -1.0, 2.0)
        self.outputGaussLegendre(sp.exp(-x ** 2 / 2), -3.0, 3.0)

if __name__ == "__main__":
    # pytest.main(["-s", "test_ch5.py::TestNumericalDifferentiation::testNumericalDifferentiation"])
    # pytest.main(["-s", "test_ch5.py::TestNewtonCotes::testNewtonCotes"])
    # pytest.main(["-s", "test_ch5.py::TestCompositeNewtonCotes::testCompositeNewtonCotes"])
    # pytest.main(["-s", "test_ch5.py::TestRomberg::testRomberg"])
    pytest.main(["-s", "test_ch5.py::TestGaussLegendre::testGaussLegendre"])
