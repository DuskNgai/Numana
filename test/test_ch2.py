import math
import os
import sys
sys.path.append(os.pardir)

import numpy as np
import pytest
import sympy as sp

import numana.chapter2 as ch2

class TestGaussJordan(object):
    def outputGaussJordan(self, A: np.ndarray, b: np.ndarray):
        x = sp.Matrix(ch2.GaussJordan(A).solve(b))
        A = sp.Matrix(A)
        b = sp.Matrix(b)
        print("Solved \033\13331m{}\033\1330m\n\033\13334m{}\033\1330m\n\033\13331m{}\033\1330m.".format(repr(A), repr(x), repr(b)))

    def testGaussJordan(self):
        A = np.array([
            [ 1,  2, -1],
            [ 2,  1, -2],
            [-3,  1,  1]
        ])
        b = np.array([3, 3, -6])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        # 1
        A = np.array([
            [ 2, -3],
            [ 5, -6],
        ])
        b = np.array([2, 8])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        A = np.array([
            [ 1,  2],
            [ 2,  3],
        ])
        b = np.array([-1, 1])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        A = np.array([
            [-1,  1],
            [ 3,  4],
        ])
        b = np.array([2, 15])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        # 2
        A = np.array([
            [ 2,  2, -1],
            [ 4,  1, -2],
            [-2,  1, -1]
        ])
        b = np.array([-2, 1, -3])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        A = np.array([
            [ 1,  2, -1],
            [ 0,  3,  1],
            [ 2, -1,  1]
        ])
        b = np.array([2, 4, 2])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        A = np.array([
            [ 2,  1, -4],
            [ 1, -1, -2],
            [-1,  3, -2]
        ])
        b = np.array([-7, -2, 6])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        # 3
        A = np.array([
            [ 3, -4,  5],
            [ 0,  3, -4],
            [ 0,  0,  5]
        ])
        b = np.array([2, -1, 5])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        A = np.array([
            [ 1, -2,  1],
            [ 0,  4, -3],
            [ 0,  0, -3]
        ])
        b = np.array([2, 1, 3])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        # 4
        A = np.array([
            [ 3, -4, -2],
            [ 6, -6,  1],
            [-3,  8,  2]
        ])
        b = np.array([3, 2, -1])[:, np.newaxis]
        self.outputGaussJordan(A, b)

        A = np.array([
            [ 2,  1, -1],
            [ 6,  2, -2],
            [ 4,  6, -3]
        ])
        b = np.array([2, 8, 5])[:, np.newaxis]
        self.outputGaussJordan(A, b)

class TestLU(object):
    def outputLUDecomposition(self, A: np.ndarray):
        lu = ch2.LU(A)
        A = sp.Matrix(A)
        L = sp.Matrix(lu.L)
        U = sp.Matrix(lu.U)
        print("Decompose \033\13331m{}\033\1330m\ninto \033\13334m{}\nand {}\033\1330m.".format(repr(A), repr(L), repr(U)))

    def outputLUSolve(self, A: np.ndarray, b: np.ndarray):
        lu = ch2.LU(A)
        x = sp.Matrix(lu.solve(b))
        L = sp.Matrix(lu.L)
        U = sp.Matrix(lu.U)
        b = sp.Matrix(b)
        print("Solved \033\13331m{}\n{}\033\1330m\n\033\13334m{}\033\1330m\n\033\13331m{}\033\1330m.".format(repr(L), repr(U), repr(x), repr(b)))

    def testLU(self):
        A = np.array([
            [2, 1, 1, 0],
            [4, 3, 3, 1],
            [8, 7, 9, 5],
            [6, 7, 9, 8]
        ])
        self.outputLUDecomposition(A)

        A = np.array([
            [-1, 4, 6],
            [-3, 14, 25],
            [1, 0, 13]
        ])
        self.outputLUDecomposition(A)

        A = np.array([[1, 2], [3, 4]])
        self.outputLUDecomposition(A)

        A = np.array([[1, 3], [2, 2]])
        self.outputLUDecomposition(A)

        A = np.array([[3, -4], [-5, 2]])
        self.outputLUDecomposition(A)

        A = np.array([
            [3, 1, 2],
            [6, 3, 4],
            [3, 1, 5]
        ])
        self.outputLUDecomposition(A)

        A = np.array([
            [4, 2, 0],
            [4, 4, 2],
            [2, 2, 3]
        ])
        self.outputLUDecomposition(A)

        A = np.array([
            [1, -1, 1, 2],
            [0, 2, 1, 0],
            [1, 3, 4, 4],
            [0, 2, 1, -1]
        ])
        self.outputLUDecomposition(A)

        A = np.array([[3, 7], [6, 1]])
        b = np.array([1, -11])
        self.outputLUSolve(A, b)

        A = np.array([[2, 3], [4, 7]])
        b = np.array([1, 3])
        self.outputLUSolve(A, b)

        A = np.array([
            [3, 1, 2],
            [6, 3, 4],
            [3, 1, 5]
        ])
        b = np.array([0, 1, 3])
        self.outputLUSolve(A, b)

        A = np.array([
            [4, 2, 0],
            [4, 4, 2],
            [2, 2, 3]
        ])
        b = np.array([2, 4, 6])
        self.outputLUSolve(A, b)

class TestCholesky(object):
    def outputCholesky(self, A: np.ndarray):
        c = ch2.Cholesky(A)
        A = sp.Matrix(A)
        R = sp.Matrix(c.R)
        print("Decompose \033\13331m{}\033\1330m\ninto \033\13334m{}\033\1330m.".format(repr(A), repr(R)))

    def testCholesky(self):
        self.A = np.array([
            [1, -1, 2, 0],
            [-1, 5, -8, 2],
            [2, -8, 14, -1],
            [0, 2, -1, 14]
        ])
        self.outputCholesky(self.A)

        self.A = np.array([
            [1, -1, 2],
            [-1, 5, 2],
            [2, 2, 17]
        ])
        self.outputCholesky(self.A)

if __name__ == "__main__":
    pytest.main(["-s", "test_ch2.py::TestGaussJordan::testGaussJordan"])
    pytest.main(["-s", "test_ch2.py::TestLU::testLU"])
    pytest.main(["-s", "test_ch2.py::TestCholesky::testCholesky"])
