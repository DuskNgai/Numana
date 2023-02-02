import math
import os
import sys
sys.path.append(os.pardir)

import pytest
import numpy as np
import sympy as sp

import numana.chapter4 as ch4

class TestLeastSquare(object):
    def setup(self):
        self.xx = sp.Symbol('x')

    def outputLeastSquare(self, A: np.ndarray, b: np.ndarray):
        x, err = ch4.LeastSquare.solve(A, b)
        for i in range(x.shape[1]):
            y = sp.Symbol(str(x[0, i]))
            for j in range(1, x.shape[0]):
                y = y + x[j, i] * self.xx ** j
            print("The \033\13331m{}\033\1330m polynomial is \033\13333m[y = {}]\033\1330m, with RMSE \033\13334m[{}]\033\1330m.".format(i + 1, y, err[i]))

    def testLeastSquare(self):
        A = np.array([[1, 1], [1, -1], [1, 1]])
        b = np.array([[2], [1], [3]])
        self.outputLeastSquare(A, b)

        A = np.array([[1, -4], [2, 3], [2, 2]])
        b = np.array([[-3], [15], [9]])
        self.outputLeastSquare(A, b)

        A = np.array([[1, -2], [1, -1], [1, 0], [1, 1], [1, 2]])
        b = np.array([[0.5], [0.5], [1.5], [3.5], [6.5]])
        self.outputLeastSquare(A, b)

        A = np.array([[1, 2], [0, -1], [-1, 0]])
        b = np.array([[-1], [3], [2]])
        self.outputLeastSquare(A, b)

        A = np.stack([np.ones(400), np.linspace(-2.0 * np.pi, 2.0 * np.pi, 400)], axis=-1)
        b = np.sin(np.linspace(-2.0 * np.pi, 2.0 * np.pi, 400).reshape(-1, 1))
        self.outputLeastSquare(A, b)


if __name__ == "__main__":
    pytest.main(["-s", "test_ch4.py::TestLeastSquare::testLeastSquare"])
