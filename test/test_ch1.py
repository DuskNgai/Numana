import os
import sys
sys.path.append(os.pardir)

import numpy as np
import pytest
import sympy as sp

import numana.chapter1 as ch1

class TestBisection(object):
    def outputBisection(self, f: sp.Function, a: float, b: float) -> None:
        solver = ch1.Solver(f)
        result = solver.bisection(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(f, a, b))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def testBisection(self) -> None:
        x = sp.Symbol('x')

        self.outputBisection(x ** 3 + x - 1, 0.0, 1.0)
        self.outputBisection(sp.cos(x) - x, 0.0, 1.0)
        # (1) (3)
        self.outputBisection(x ** 3 - 9, 2.0, 3.0)
        self.outputBisection(3 * x ** 3 + x ** 2 - x - 5 - 9, 1.0, 2.0)
        self.outputBisection(sp.cos(x) ** 2 - x + 6, 6.0, 7.0)
        # (2) (4)
        self.outputBisection(x ** 5 + x - 1, 0.0, 1.0)
        self.outputBisection(sp.sin(x) - 6 * x - 5, -1.0, 0.0)
        self.outputBisection(sp.log(x) + x ** 2 - 3, 1.0, 2.0)
        # (5)
        self.outputBisection(x ** 4 - x ** 3 - 10, 2.0, 3.0)

        # (3)
        self.outputBisection(2 * x ** 3 - 6 * x - 1, -1.0, 0.0)
        self.outputBisection(sp.exp(x - 2) + x ** 3 - x, -0.5, 0.5)
        self.outputBisection(1 + 5 * x - 6 * x ** 3 + sp.exp(x * 2), -0.5, 0.5)
        # (6)
        self.outputBisection(sp.cos(x) - sp.sin(x), 0.0, 1.0)
        # (7)
        self.outputBisection(
            sp.Matrix([[1, 2, 3, x], [4, 5, x, 6], [7, x, 8, 9], [x, 10, 11, 12]]).det() - 1000,
            -18.0, -17.0
        )
        self.outputBisection(
            sp.Matrix([[1, 2, 3, x], [4, 5, x, 6], [7, x, 8, 9], [x, 10, 11, 12]]).det() - 1000,
            9.0, 10.0
        )
        # (8)
        # (9)
        self.outputBisection(sp.pi * x ** 2 * (1 - 1 / 3 * x) - 1, 0.0, 1.0)

class TestFixedPointIteration(object):
    def outputFixedPointIteration(self, f: sp.Function, a: float) -> None:
        solver = ch1.Solver(f)
        result = solver.fixedPointIteration(a)
        print("The function is \033\13331mf(x) = {}\033\1330m and the initial value is \033\13331m[{}]\033\1330m.".format(f, a))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def testFixedPointIteration(self) -> None:
        x = sp.Symbol('x')
        self.outputFixedPointIteration(1 - x ** 3, 0.5)
        self.outputFixedPointIteration((1 - x) ** (1 / 3), 0.5)
        self.outputFixedPointIteration((1 + 2 * x ** 3) / (1 + 3 * x ** 2), 0.5)
        self.outputFixedPointIteration((2 / 3) * x + (8 / 3)/ (x ** 2), 4.0)
        # (1)
        self.outputFixedPointIteration((x ** 2 + 2 * x + 2) / (x ** 2 + x), 1.0)
        self.outputFixedPointIteration(sp.log(7 - x), 1.0)
        self.outputFixedPointIteration(sp.log(4 - sp.sin(x)), 1.0)
        # (2)
        self.outputFixedPointIteration((x ** 4 + x ** 3 + x ** 2 + 1) / (x ** 4 + x ** 3 + x ** 2 + x + 1), 0.0)
        self.outputFixedPointIteration((sp.sin(x) - 5) / 6, 0.0)
        self.outputFixedPointIteration(sp.sqrt(3 - sp.log(x)), 2.0)
        # (3)
        self.outputFixedPointIteration((x + 3 / x) / 2, 2.0)
        self.outputFixedPointIteration((x + 5 / x) / 2, 2.0)
        # (4)
        self.outputFixedPointIteration((2 / 3) * x + (2 / 3) / (x ** 2), 1.0)
        self.outputFixedPointIteration((2 / 3) * x + (3 / 3) / (x ** 2), 1.0)
        self.outputFixedPointIteration((2 / 3) * x + (5 / 3) / (x ** 2), 1.0)
        # (5)
        self.outputFixedPointIteration((sp.cos(x) ** 2 + x ** 2) / (x + 1), 0.1)
        # (7)
        self.outputFixedPointIteration(1 - 5 * x + 15 / 2 * x ** 2 - 5 / 2 * x ** 3, 2.18)

class TestNewton(object):
    def outputNewton(self, f: sp.Function, a: float) -> None:
        solver = ch1.Solver(f)
        result = solver.newton(a)
        print("The function is \033\13331mf(x) = {}\033\1330m and the initial value is \033\13331m[{}]\033\1330m.".format(f, a))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def testNewton(self):
        x = sp.Symbol('x')
        # (1)
        self.outputNewton(x ** 3 + x - 2, 0.0)
        self.outputNewton(x ** 4 - x ** 2 + x - 1, 0.0)
        self.outputNewton(x ** 2 - x - 1, 0.0)
        # (2)
        self.outputNewton(x ** 3 + x ** 2 - 1, 1.0)
        self.outputNewton(x ** 2 + 1 / (x + 1) - 3 * x, 1.0)
        self.outputNewton(5 * x - 10, 1.0)
        # (12) inf
        # self.outputNewton(1 / x, 1.0)
        # (13)
        self.outputNewton(x ** 3 - 4 * x, 0.5)
        # (1)
        self.outputNewton(x ** 3 - 2 * x - 2, 1.0)
        self.outputNewton(sp.exp(x) + x - 7, 1.0)
        self.outputNewton(sp.exp(x) + sp.sin(x) - 4, 1.0)
        # (2)
        self.outputNewton(x ** 5 + x - 1, 1.0)
        self.outputNewton(sp.sin(x) - 6 * x - 5, 1.0)
        self.outputNewton(sp.log(x) + x ** 2 - 3, 1.0)
        # (3)
        self.outputNewton(27 * x ** 3 + 54 * x ** 2 + 36 * x + 8, -1.0)

class TestSecant(object):
    def outputSecant(self, f: sp.Function, a: float, b: float) -> None:
        solver = ch1.Solver(f)
        result = solver.secant(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(f, a, b))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def testSecant(self) -> None:
        x = sp.Symbol('x')
        self.outputSecant(x ** 3 + x - 1, 0.0, 1.0)
        self.outputSecant(x ** 3 - 2 * x - 2, 1.0, 2.0)
        self.outputSecant(sp.exp(x) + x - 7, 1.0, 2.0)
        self.outputSecant(sp.exp(x) + sp.sin(x) - 4, 1.0, 2.0)

class TestRegulaFalsi(object):
    def outputRegulaFalsi(self,f : sp.Function, a: float, b: float) -> None:
        solver = ch1.Solver(f)
        result = solver.regulaFalsi(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(f, a, b))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def testRegulaFalsi(self) -> None:
        x = sp.Symbol('x')
        self.outputRegulaFalsi(x ** 3 - 2 * x - 2, 1.0, 2.0)
        self.outputRegulaFalsi(sp.exp(x) + x - 7, 1.0, 2.0)
        self.outputRegulaFalsi(sp.exp(x) + sp.sin(x) - 4, 1.0, 2.0)

class TestInverseInterpolation(object):
    def outputInverseInterpolation(self, f: sp.Function, a: float, b: float, c: float) -> None:
        solver = ch1.Solver(f)
        result = solver.inverseInterpolation(a, b, c)
        print("The function is \033\13331mf(x) = {}\033\1330m and the initials are \033\13331m[{}, {}, {}]\033\1330m.".format(f, a, b, c))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def testInverseInterpolation(self) -> None:
        x = sp.Symbol('x')
        self.outputInverseInterpolation(x ** 3 - 2 * x - 2, 1.0, 2.0, 0.0)
        self.outputInverseInterpolation(sp.exp(x) + x - 7, 1.0, 2.0, 0.0)
        self.outputInverseInterpolation(sp.exp(x) + sp.sin(x) - 4, 1.0, 2.0, 0.0)

if __name__ == "__main__":
    pytest.main(["-s", "test_ch1.py::TestBisection::testBisection"])
    pytest.main(["-s", "test_ch1.py::TestFixedPointIteration::testFixedPointIteration"])
    pytest.main(["-s", "test_ch1.py::TestNewton::testNewton"])
    pytest.main(["-s", "test_ch1.py::TestSecant::testSecant"])
    pytest.main(["-s", "test_ch1.py::TestRegulaFalsi::testRegulaFalsi"])
    pytest.main(["-s", "test_ch1.py::TestInverseInterpolation::testInverseInterpolation"])
