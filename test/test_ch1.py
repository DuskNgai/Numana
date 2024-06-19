from pathlib import Path
import sys

import pytest
import sympy as sp

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numana.chapter1 as ch1

class TestBisection(object):
    def output_bisection(self, f: sp.Function, a: float, b: float) -> None:
        solver = ch1.Solver(f)
        result = solver.bisection(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(f, a, b))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def test_bisection(self) -> None:
        x = sp.symbols('x')

        self.output_bisection(x ** 3 + x - 1, 0.0, 1.0)
        self.output_bisection(sp.cos(x) - x, 0.0, 1.0)
        # (1) (3)
        self.output_bisection(x ** 3 - 9, 2.0, 3.0)
        self.output_bisection(3 * x ** 3 + x ** 2 - x - 5 - 9, 1.0, 2.0)
        self.output_bisection(sp.cos(x) ** 2 - x + 6, 6.0, 7.0)
        # (2) (4)
        self.output_bisection(x ** 5 + x - 1, 0.0, 1.0)
        self.output_bisection(sp.sin(x) - 6 * x - 5, -1.0, 0.0)
        self.output_bisection(sp.log(x) + x ** 2 - 3, 1.0, 2.0)
        # (5)
        self.output_bisection(x ** 4 - x ** 3 - 10, 2.0, 3.0)

        # (3)
        self.output_bisection(2 * x ** 3 - 6 * x - 1, -1.0, 0.0)
        self.output_bisection(sp.exp(x - 2) + x ** 3 - x, -0.5, 0.5)
        self.output_bisection(1 + 5 * x - 6 * x ** 3 + sp.exp(x * 2), -0.5, 0.5)
        # (6)
        self.output_bisection(sp.cos(x) - sp.sin(x), 0.0, 1.0)
        # (7)
        self.output_bisection(
            sp.Matrix([[1, 2, 3, x], [4, 5, x, 6], [7, x, 8, 9], [x, 10, 11, 12]]).det() - 1000,
            -18.0, -17.0
        )
        self.output_bisection(
            sp.Matrix([[1, 2, 3, x], [4, 5, x, 6], [7, x, 8, 9], [x, 10, 11, 12]]).det() - 1000,
            9.0, 10.0
        )
        # (8)
        # (9)
        self.output_bisection(sp.pi * x ** 2 * (1 - 1 / 3 * x) - 1, 0.0, 1.0)

class TestFixedPointIteration(object):
    def output_fixed_point_iteration(self, f: sp.Function, a: float) -> None:
        solver = ch1.Solver(f)
        result = solver.fixed_point_iteration(a)
        print("The function is \033\13331mf(x) = {}\033\1330m and the initial value is \033\13331m[{}]\033\1330m.".format(f, a))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def test_fixed_point_iteration(self) -> None:
        x = sp.symbols('x')

        self.output_fixed_point_iteration(1 - x ** 3, 0.5)
        self.output_fixed_point_iteration((1 - x) ** (1 / 3), 0.5)
        self.output_fixed_point_iteration((1 + 2 * x ** 3) / (1 + 3 * x ** 2), 0.5)
        self.output_fixed_point_iteration((2 / 3) * x + (8 / 3)/ (x ** 2), 4.0)
        # (1)
        self.output_fixed_point_iteration((x ** 2 + 2 * x + 2) / (x ** 2 + x), 1.0)
        self.output_fixed_point_iteration(sp.log(7 - x), 1.0)
        self.output_fixed_point_iteration(sp.log(4 - sp.sin(x)), 1.0)
        # (2)
        self.output_fixed_point_iteration((x ** 4 + x ** 3 + x ** 2 + 1) / (x ** 4 + x ** 3 + x ** 2 + x + 1), 0.0)
        self.output_fixed_point_iteration((sp.sin(x) - 5) / 6, 0.0)
        self.output_fixed_point_iteration(sp.sqrt(3 - sp.log(x)), 2.0)
        # (3)
        self.output_fixed_point_iteration((x + 3 / x) / 2, 2.0)
        self.output_fixed_point_iteration((x + 5 / x) / 2, 2.0)
        # (4)
        self.output_fixed_point_iteration((2 / 3) * x + (2 / 3) / (x ** 2), 1.0)
        self.output_fixed_point_iteration((2 / 3) * x + (3 / 3) / (x ** 2), 1.0)
        self.output_fixed_point_iteration((2 / 3) * x + (5 / 3) / (x ** 2), 1.0)
        # (5)
        self.output_fixed_point_iteration((sp.cos(x) ** 2 + x ** 2) / (x + 1), 0.1)
        # (7)
        self.output_fixed_point_iteration(1 - 5 * x + 15 / 2 * x ** 2 - 5 / 2 * x ** 3, 2.18)

class TestNewton(object):
    def output_newton(self, f: sp.Function, a: float) -> None:
        solver = ch1.Solver(f)
        result = solver.newton(a)
        print("The function is \033\13331mf(x) = {}\033\1330m and the initial value is \033\13331m[{}]\033\1330m.".format(f, a))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def test_newton(self):
        x = sp.symbols('x')

        # (1)
        self.output_newton(x ** 3 + x - 2, 0.0)
        self.output_newton(x ** 4 - x ** 2 + x - 1, 0.0)
        self.output_newton(x ** 2 - x - 1, 0.0)
        # (2)
        self.output_newton(x ** 3 + x ** 2 - 1, 1.0)
        self.output_newton(x ** 2 + 1 / (x + 1) - 3 * x, 1.0)
        self.output_newton(5 * x - 10, 1.0)
        # (12) inf
        # self.output_newton(1 / x, 1.0)
        # (13)
        self.output_newton(x ** 3 - 4 * x, 0.5)
        # (1)
        self.output_newton(x ** 3 - 2 * x - 2, 1.0)
        self.output_newton(sp.exp(x) + x - 7, 1.0)
        self.output_newton(sp.exp(x) + sp.sin(x) - 4, 1.0)
        # (2)
        self.output_newton(x ** 5 + x - 1, 1.0)
        self.output_newton(sp.sin(x) - 6 * x - 5, 1.0)
        self.output_newton(sp.log(x) + x ** 2 - 3, 1.0)
        # (3)
        self.output_newton(27 * x ** 3 + 54 * x ** 2 + 36 * x + 8, -1.0)

class TestSecant(object):
    def output_secant(self, f: sp.Function, a: float, b: float) -> None:
        solver = ch1.Solver(f)
        result = solver.secant(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(f, a, b))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def test_secant(self) -> None:
        x = sp.symbols('x')

        self.output_secant(x ** 3 + x - 1, 0.0, 1.0)
        self.output_secant(x ** 3 - 2 * x - 2, 1.0, 2.0)
        self.output_secant(sp.exp(x) + x - 7, 1.0, 2.0)
        self.output_secant(sp.exp(x) + sp.sin(x) - 4, 1.0, 2.0)

class TestRegulaFalsi(object):
    def output_regula_falsi(self,f : sp.Function, a: float, b: float) -> None:
        solver = ch1.Solver(f)
        result = solver.regula_falsi(a, b)
        print("The function is \033\13331mf(x) = {}\033\1330m and the interval is \033\13331m[{}, {}]\033\1330m.".format(f, a, b))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def test_regula_falsi(self) -> None:
        x = sp.symbols('x')

        self.output_regula_falsi(x ** 3 - 2 * x - 2, 1.0, 2.0)
        self.output_regula_falsi(sp.exp(x) + x - 7, 1.0, 2.0)
        self.output_regula_falsi(sp.exp(x) + sp.sin(x) - 4, 1.0, 2.0)

class TestInverseInterpolation(object):
    def output_inverse_interpolation(self, f: sp.Function, a: float, b: float, c: float) -> None:
        solver = ch1.Solver(f)
        result = solver.inverse_interpolation(a, b, c)
        print("The function is \033\13331mf(x) = {}\033\1330m and the initials are \033\13331m[{}, {}, {}]\033\1330m.".format(f, a, b, c))
        print("The result is: \033\13334m[{}]\033\1330m.".format(result))

    def test_inverse_interpolation(self) -> None:
        x = sp.symbols('x')

        self.output_inverse_interpolation(x ** 3 - 2 * x - 2, 1.0, 2.0, 0.0)
        self.output_inverse_interpolation(sp.exp(x) + x - 7, 1.0, 2.0, 0.0)
        self.output_inverse_interpolation(sp.exp(x) + sp.sin(x) - 4, 1.0, 2.0, 0.0)

if __name__ == "__main__":
    pytest.main(["-s", "test_ch1.py::TestBisection::test_bisection"])
    pytest.main(["-s", "test_ch1.py::TestFixedPointIteration::test_fixed_point_iteration"])
    pytest.main(["-s", "test_ch1.py::TestNewton::test_newton"])
    pytest.main(["-s", "test_ch1.py::TestSecant::test_secant"])
    pytest.main(["-s", "test_ch1.py::TestRegulaFalsi::test_regula_falsi"])
    pytest.main(["-s", "test_ch1.py::TestInverseInterpolation::test_inverse_interpolation"])
