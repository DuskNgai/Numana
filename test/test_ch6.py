import math
import os
import sys
sys.path.append(os.pardir)

import pytest
import sympy as sp
import sympy.abc

import numana.chapter6 as ch6

class TestExplicitEuler(object):
    def outputExplicitEuler(self, f: sp.Function, initial: float, num_steps: int, start: float, end: float):
        solver = ch6.EulerMethod(f)
        result = solver.explicit(initial, num_steps, (start, end))
        print(result)

    def testExplicitEuler(self):
        x = sympy.abc.x
        y = sympy.abc.y

        self.outputExplicitEuler(x * y + x ** 3, 1.0, 5, 0.0, 1.0)

    def outputTrapezoidEuler(self, f: sp.Function, initial: float, num_steps: int, start: float, end: float):
        solver = ch6.EulerMethod(f)
        result = solver.trapezoid(initial, num_steps, (start, end))
        print(result)

    def testTrapezoidEuler(self):
        x = sympy.abc.x
        y = sympy.abc.y

        self.outputTrapezoidEuler(x * y + x ** 3, 1.0, 10, 0.0, 1.0)

if __name__ == "__main__":
    # pytest.main(["-s", "test_ch6.py::TestExplicitEuler::testExplicitEuler"])
    pytest.main(["-s", "test_ch6.py::TestExplicitEuler::testTrapezoidEuler"])
