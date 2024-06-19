from pathlib import Path
import sys

import pytest
import numpy as np
import sympy as sp
import sympy.abc

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numana.chapter4 as ch4

class TestLeastSquare(object):
    def setup_method(self):
        self.x = sympy.abc.x

    def output_least_square(self, A: np.ndarray, b: np.ndarray):
        x, err = ch4.LeastSquare.solve(A, b)
        for i in range(x.shape[1]):
            y = sp.Symbol(str(x[0, i]))
            for j in range(1, x.shape[0]):
                y = y + x[j, i] * self.x ** j
            print("The \033\13331m{}\033\1330m polynomial is \033\13333m[y = {}]\033\1330m, with RMSE \033\13334m[{}]\033\1330m.".format(i + 1, y, err[i]))

    def test_least_square(self):
        A = np.array([[1, 1], [1, -1], [1, 1]])
        b = np.array([[2], [1], [3]])
        self.output_least_square(A, b)

        A = np.array([[1, -4], [2, 3], [2, 2]])
        b = np.array([[-3], [15], [9]])
        self.output_least_square(A, b)

        A = np.array([[1, -2], [1, -1], [1, 0], [1, 1], [1, 2]])
        b = np.array([[0.5], [0.5], [1.5], [3.5], [6.5]])
        self.output_least_square(A, b)

        A = np.array([[1, 2], [0, -1], [-1, 0]])
        b = np.array([[-1], [3], [2]])
        self.output_least_square(A, b)

        A = np.stack([np.ones(400), np.linspace(-2.0 * np.pi, 2.0 * np.pi, 400)], axis=-1)
        b = np.sin(np.linspace(-2.0 * np.pi, 2.0 * np.pi, 400).reshape(-1, 1))
        self.output_least_square(A, b)


if __name__ == "__main__":
    pytest.main(["-s", "test_ch4.py::TestLeastSquare::test_least_square"])
