from pathlib import Path
import sys

import pytest
import sympy as sp

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numana.chapter6 as ch6

class TestExplicitEuler(object):
    def output_explicit_euler(self, f: sp.Function, initial: float, num_steps: int, start: float, end: float):
        solver = ch6.EulerMethod(f)
        result = solver.explicit(initial, num_steps, (start, end))
        print(result)

    def test_explicit_euler(self):
        x = sp.symbols('x')
        y = sp.symbols('y')

        self.output_explicit_euler(x * y + x ** 3, 1.0, 5, 0.0, 1.0)

    def output_trapezoid_euler(self, f: sp.Function, initial: float, num_steps: int, start: float, end: float):
        solver = ch6.EulerMethod(f)
        result = solver.trapezoid(initial, num_steps, (start, end))
        print(result)

    def test_trapezoid_euler(self):
        x = sp.symbols('x')
        y = sp.symbols('y')

        self.output_trapezoid_euler(x * y + x ** 3, 1.0, 10, 0.0, 1.0)

if __name__ == "__main__":
    pytest.main(["-s", "test_ch6.py::TestExplicitEuler::test_explicit_euler"])
    pytest.main(["-s", "test_ch6.py::TestExplicitEuler::test_trapezoid_euler"])
