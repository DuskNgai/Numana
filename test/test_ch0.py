import math
import os
import sys
from typing import Optional, Union
sys.path.append(os.pardir)

import numpy as np
import pytest
import sympy as sp
import sympy.abc

import numana.chapter0 as ch0

class TestNest(object):
    def setup(self):
        self.x = sympy.abc.x

    def outputNest(self,
        x: Union[float, list[float]],
        coefficients: list[float],
        base_points: Optional[list[float]] = None
    ):
        coeffs = np.asarray(coefficients[-2::-1])
        if base_points is None:
            base_points = np.zeros_like(coeffs)
        else:
            base_points = np.asarray(base_points)

        y = ch0.nest(x, coefficients, base_points)
        yy = sp.Symbol(str(coefficients[-1]))

        for i in range(coeffs.shape[0]):
            yy = yy * (self.x + base_points[i]) + coeffs[i]
        yy = sp.expand(yy, mul=True)
        print("Evaluating \033\13331m[{}]\033\1330m at x = {}, the value is \033\13334m[{}]\033\1330m.".format(yy, x, y))
        
    def testNest(self):
        self.outputNest(0.5,   [-1, 5, -3, 3, 2],   [0, 0, 0, 0])
        self.outputNest(0.5,   [-1, 5, -3, 3, 2])
        self.outputNest([-2, -1, 0, 1, 2], [-1, 5, -3, 3, 2])
        self.outputNest(1,     [1, 0.5, 0.5, -0.5], [0, 2, 3])

        self.outputNest(1 / 3, [1, 1, 5, 1, 6])
        self.outputNest(1 / 3, [1, -5, 5, 4, -3])
        self.outputNest(1 / 3, [1, 0, -1, 1, 2])

        self.outputNest(-0.5,  [7, -3, -2, 6])
        self.outputNest(-0.5,  [1, -3, 1, -3, -1, 8])
        self.outputNest(-0.5,  [4, -2, 0, 0, -2, 0, 4])

        self.outputNest(0.5,   [1, 0, 2, 0, -4, 0, 1])
        self.outputNest(0.25,  [1, 2, -4, 1])

        self.outputNest(5,     [1, 0.5, 0.5, -0.5], [0, 2, 3])
        self.outputNest(-1,    [1, 0.5, 0.5, -0.5], [0, 2, 3])
        self.outputNest(0.5,   [4, 4, 1, 3, 2],     [0, 1, 2, 3])
        self.outputNest(-0.5,  [4, 4, 1, 3, 2],     [0, 1, 2, 3])

        self.outputNest(1.00001, np.ones(51, dtype=int))
        print("\033\13334m[{}]\033\1330m.".format((1.00001 ** 51 - 1) / (1.00001 - 1)))

        self.outputNest(1.00001, np.array([1, -1] * 50))
        print("\033\13334m[{}]\033\1330m.".format((1 - 1.00001) * ch0.nest(1.00001 ** 2, np.ones(50, dtype=int))))

class TestSignificance(object):
    def setup(self):
        self.x = sympy.abc.x

    def outputQuadratic(self, a: float, b: float, c: float):
        solution = ch0.quadratic(a, b, c)
        print("Solving \033\13331m[{}]\033\1330m in real number, the solution is \033\13334m{}\033\1330m.".format(a * self.x ** 2 + b * self.x + c, solution))

    def testQuadratic(self):
        self.outputQuadratic(1, -3, 2)
        self.outputQuadratic(1, 3, 8 ** -14)
        self.outputQuadratic(1, 100, -1e-12)

    def testSignificance(self):
        x = np.logspace(-1, -14, 14)
        print((1 - 1 / np.cos(x)) / (np.tan(x) * np.tan(x)))
        print(-(1 / (1 + 1 / np.cos(x))))
        print((1 - (1 - x) ** 3) / x)
        print(x ** 2 - 3 * x + 3)

        # taylor expansion
        print((np.tan(x) - x) / (x ** 3))
        print(1 / 3 + 2 / 15 * x ** 2 + 17 / 315 * x ** 4)
        print((np.exp(x) + np.cos(x) - np.sin(x) - 2) / (x ** 3))
        print(1 / 3 + x / 12 + x ** 4 / 2520)

        print(math.sqrt(12345678987654321 ** 2 + 123 ** 2) - 12345678987654321)
        print((123 ** 2) / (math.sqrt(12345678987654321 ** 2 + 123 ** 2) + 12345678987654321))
        print(math.sqrt(246886422468 ** 2 + 13579) - 246886422468)
        print(13579 / (math.sqrt(246886422468 ** 2 + 13579) + 246886422468))
        print(math.sqrt(3344556600 ** 2 + 1.2222222 ** 2) - 3344556600)
        print((1.2222222 ** 2) / (math.sqrt(3344556600 ** 2 + 1.2222222 ** 2) + 3344556600))

class TestBinaryDecimal(object):
    def outputDec2bin(self, deci: Union[float, str], fraction_len: int = 10):
        result = ch0.dec2bin(deci, fraction_len)
        print("Convert \033\13331m[{}]\033\1330m into \033\13334m[{}]\033\1330m.".format(deci, result))

    def testDec2bin(self):
        self.outputDec2bin(64)
        self.outputDec2bin(17)
        self.outputDec2bin(79)
        self.outputDec2bin(227)

        self.outputDec2bin(1 / 8)
        self.outputDec2bin(7 / 8)
        self.outputDec2bin(35 / 16)
        self.outputDec2bin(31 / 64)

        self.outputDec2bin(10.5)
        self.outputDec2bin(1 / 3)
        self.outputDec2bin(5 / 7)
        self.outputDec2bin(12.8)
        self.outputDec2bin(55.4)
        self.outputDec2bin(0.1)

        self.outputDec2bin(11.25)
        self.outputDec2bin(2 / 3)
        self.outputDec2bin(3/ 5)
        self.outputDec2bin(3.2)
        self.outputDec2bin(30.6)
        self.outputDec2bin(99.9)

        self.outputDec2bin(3.1415926535897932, 15)
        self.outputDec2bin(2.7182818284590452, 15)

        self.outputDec2bin(1 / 4, 53)
        self.outputDec2bin(1 / 3, 53)
        self.outputDec2bin(2 / 3, 53)
        self.outputDec2bin(0.9, 53)
        self.outputDec2bin(9.6, 53)
        self.outputDec2bin(100.2, 53)
        self.outputDec2bin(44 / 7, 53)
        self.outputDec2bin(7 / 3, 53)
        self.outputDec2bin(4 / 3, 53)
        self.outputDec2bin(3.3, 53)
        self.outputDec2bin(9 / 7, 53)
        self.outputDec2bin(2.75, 53)
        self.outputDec2bin(2.7, 53)
        self.outputDec2bin(10 / 3, 53)

    def outputBin2dec(self, bina: str, recur_pos: tuple = None):
        result = ch0.bin2dec(bina, recur_pos)
        print("Convert \033\13331m[{}]\033\1330m into \033\13334m[{}]\033\1330m.".format(bina, result))

    def testBin2dec(self):
        self.outputBin2dec('1010101')
        self.outputBin2dec('1011.101')
        self.outputBin2dec('10111.01', (0, 1))
        self.outputBin2dec('110.10', (0, 1))
        self.outputBin2dec('10.110', (0, 2))
        self.outputBin2dec('110.1101', (1, 3))
        self.outputBin2dec('10.0101101', (3, 6))
        self.outputBin2dec('111.1', (0, 0))

        self.outputBin2dec('11011')
        self.outputBin2dec('110111.001')
        self.outputBin2dec('111.001', (0, 2))
        self.outputBin2dec('1010.01', (0, 1))
        self.outputBin2dec('10111.10', (0, 1))
        self.outputBin2dec('1111.010001', (3, 5))

if __name__ == "__main__":
    pytest.main(["-s", "test_ch0.py::TestNest::testNest"])
    pytest.main(["-s", "test_ch0.py::TestSignificance::testQuadratic"])
    pytest.main(["-s", "test_ch0.py::TestSignificance::testSignificance"])
    pytest.main(["-s", "test_ch0.py::TestBinaryDecimal::testDec2bin"])
    pytest.main(["-s", "test_ch0.py::TestBinaryDecimal::testBin2dec"])
