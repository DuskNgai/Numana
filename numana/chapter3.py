from numbers import Number

import numpy as np
import matplotlib.pyplot as plt


class LagrangeInterpolation:
    """
    Interpolates input points, work as a functor.
    @param x, y: the points to interpolate.
    """
    def __init__(self, x: np.ndarray, y: np.ndarray):
        self.x = np.asarray(x, dtype=float).reshape(1, -1)
        self.y = np.asarray(y, dtype=float).reshape(-1)

        assert self.x.shape[-1] == self.y.shape[-1], "x and y must be the length"
        assert np.all(self.x[:-1] < self.x[1:]), "x should be strictly increasing"

        self.n = self.x.shape[-1]

        self.xi_xj = self.x.T - self.x
        print(self.xi_xj)

    def __call__(self, x0: Number | np.ndarray) -> Number | np.ndarray:
        if isinstance(x0, (int, float)):
            x0 = float(x0)
            # numpy version, only for single x0
            x0_xj = x0 - np.tile(self.x, (self.n, 1))
            # let diagonal element be 1 after divided by x0_xj
            np.fill_diagonal(self.xi_xj, np.diag(x0_xj))
            w = x0_xj / self.xi_xj
            return np.dot(np.prod(w, axis=1), self.y)
        else:
            x0 = np.array(x0, dtype=float).reshape(-1)
            # python version, also use numpy for multiple x0
            res = 0.0
            for i in range(self.n):
                w = 1.0
                for j in range(self.n):
                    if i != j:
                        w *= (x0 - self.x[0, j]) / self.xi_xj[i, j]
                res += w * self.y[i]
            return res

class NewtonDividedDifference:
    """
    Using Newton's divided difference to interpolate points, work as a functor.
    @param x, y: the points to interpolate.
    """

    def __init__(self, x: np.ndarray, y: np.ndarray):
        self.x = np.asarray(x, dtype=float).reshape(1, -1)
        self.y = np.asarray(y, dtype=float).reshape(1, -1)

        assert self.x.shape[-1] == self.y.shape[-1], "x and y must be the length"
        assert np.all(self.x[:-1] < self.x[1:]), "x should be strictly increasing"

        self.n = self.x.shape[-1]

        self.xi_xj = np.zeros((self.n, self.n))
        for i in range(0, self.n):
            self.xi_xj[:self.n - i, i] = self.x[0, i:]
        self.xi_xj = self.xi_xj - self.x.T

        # arrange all of the coefficients in upper triangular form
        self.coefficients = np.zeros((self.n, self.n))
        self.coefficients[:, 0] = self.y
        for i in range(1, self.n):
            self.coefficients[:self.n - i, i] = np.diff(self.coefficients[:, i - 1])[:self.n - i] / self.xi_xj[:self.n - i, i]

    def __call__(self, x0: Number | np.ndarray) -> Number | np.ndarray:
        if isinstance(x0, (int, float)):
            x0 = float(x0)
        else:
            x0 = np.array(x0, dtype=float).reshape(-1)
        # numpy version evaluation
        res = 0.0
        for i in range(self.n - 1, -1, -1):
            res = res * (x0 - self.x[0, i]) + self.coefficients[0, i]
        return res

class CubicSpline:
    """
    Using Cubic Splines to fit the given point, often used in plotting, work as a functor.
    @param x, y: the points to interpolate.
    @param edge_case: natural, clamp, not-a-knot
    """
    def __init__(self, x: np.ndarray, y: np.ndarray, edge_case: str = 'natural'):
        self.x = np.asarray(x, dtype=float).reshape(-1)
        self.y = np.asarray(y, dtype=float).reshape(-1)

        assert self.x.shape[-1] == self.y.shape[-1], "x and y must be the length"
        assert np.all(self.x[:-1] < self.x[1:]), "x should be strictly increasing"

        self.n = self.x.shape[0]

        dx = self.x[1:] - self.x[:-1]
        dy = self.y[1:] - self.y[:-1]

        A = np.zeros((self.n, self.n), dtype=float)
        for i in range(1, self.n - 1):
            A[i, i - 1:i + 2] = np.array([dx[i - 1], 2 * (dx[i - 1] + dx[i]), dx[i]], dtype=float)

        B = np.zeros(self.n)
        B[1:-1] = 3 * ((dy[1:] / dx[1:]) - (dy[:-1] / dx[:-1]))

        if edge_case == 'natural':
            assert self.n >= 2, "natural edge case needs at least 2 points"
            A[0, 0] = 1.0
            A[-1, -1] = 1.0
        elif edge_case == 'clamp':
            assert self.n >= 3, "clamp edge case needs at least 3 points"
            A[0, 0:2] = np.array([1.0, -1.0], dtype=float)
            A[-1, -2:] = np.array([1.0, -1.0], dtype=float)
        elif edge_case == 'not-a-knot':
            assert self.n >= 4, "not-a-knot edge case needs at least 4 points"
            A[0, 0:3] = np.array([-dx[1], dx[1] + dx[0], -dx[0]], dtype=float)
            A[-1, -3:] = np.array([-dx[-1], dx[-1] + dx[-2], -dx[-2]], dtype=float)
        else:
            raise TypeError("Not a valid edge case.")

        self.c = np.linalg.solve(A, B)
        self.d = (self.c[1:] - self.c[:-1]) / (3 * dx)
        self.b = (dy / dx) - (dx / 3) * (2 * self.c[:-1] + self.c[1:])

    def __call__(self, x0: Number | np.ndarray) -> Number | np.ndarray:
        if isinstance(x0, (int, float)):
            x0 = float(x0)
        else:
            x0 = np.array(x0, dtype=float).reshape(-1)

        # binary search
        assert np.min(x0) >= self.x[0] and np.max(x0) <= self.x[-1], "Invalid evaluation"
        index = np.searchsorted(self.x, x0, side='left') - 1
        # fancy indexing only for array-like objects
        if isinstance(x0, np.ndarray):
            index[index == -1] = 0
        index = index.astype(np.int32)

        x0_xi = x0 - self.x[index]
        return (((self.d[index] * x0_xi) + self.c[index]) * x0_xi + self.b[index]) * x0_xi + self.y[index]

class Bezier:
    """
    Using Bezier curve to fit the given point, often used in plotting, work as a functor.
    @param x, y: the points to interpolate.
    Pass an array with the range [0, 1] when call it as a function.
    """
    def __init__(self, x: Iterable, y: Iterable):
        self.x = np.asarray(x, dtype=float).reshape(-1)
        self.y = np.asarray(y, dtype=float).reshape(-1)

        assert self.x.shape[0] == self.y.shape[0], "x and y must be the length"

    def __call__(self, t: Union[float, Iterable]):
        self.t = np.asarray(t, dtype=float).reshape(-1)
        # recursive maybe faster
        return self.__bezierCurve(self.x, self.y)

    def __bezierCurve(self, xs: np.ndarray, ys: np.ndarray):
        if xs.shape[0] == 1:
            return xs, ys
        else:
            return (1 - self.t) * self.__bezierCurve(xs[:-1], ys[:-1]) + self.t * self.__bezierCurve(xs[1:], ys[1:])

class BSpline(Bezier):
    def __init__(self, x: Iterable, y: Iterable):
        super(BSpline, self).__init__(x, y)

    def __call__(self, t: Union[float, Iterable]):
        self.t = np.asarray(t, dtype=float).reshape(-1)
        # recursive maybe faster
        return self.__bezierCurve(self.x, self.y)

if __name__ == "__main__":
    x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    y = [0, 1, -2, 2, 3, 1, -1, 2, 1, 3]

    lagrange = True

    if lagrange:
        x = [0,0.5,1]
        y = [1,np.exp(-0.5),np.exp(-1)]
        f = LagrangeInterpolation(x, y)
        ic(f(0.5))
        plt.figure(figsize=(8, 6))
        t = np.linspace(0, 1, 500)
        plt.plot(t, f(t) - np.exp(-t))
        # plt.scatter(x, y)
        plt.show()

        # g = NewtonDividedDifference(x, y)
        # ic(g(0.5))
        # plt.figure(figsize=(8, 6))
        # t = np.linspace(0, 9, 500)
        # plt.plot(t, g(t))
        # plt.scatter(x, y)
        # plt.show()

        # 17
        # year = [1800, 1850, 1900, 2000]
        # CO2 = [280, 283, 291, 370]
        # order3 = LagrangeInterpolation(year, CO2)
        # ic(order3(1950))
        # ic(order3(2050))

        # 18
        # temperature = [25, 40, 50, 60]
        # time = [95, 75, 63, 54]
        # order2 = LagrangeInterpolation(temperature[1:], time[1:])
        # order3 = LagrangeInterpolation(temperature, time)
        # ic(order2(70))
        # ic(order3(70))

        # 1
        # year = [1960, 1970, 1990, 2000]
        # population = [3039585530, 3707475887, 5281653820, 6079603571]
        # order1 = LagrangeInterpolation(year[0:2], population[0:2])
        # order2 = LagrangeInterpolation(year[0:3], population[0:3])
        # order3 = LagrangeInterpolation(year, population)
        # ic(order1(1980))
        # ic(order2(1980))
        # ic(order3(1980))

        # 4
        # u = np.array([0, np.pi / 6, np.pi / 3, np.pi / 2])
        # sinu = np.sin(u)
        # cosu = np.cos(u)
        # my_sin = NewtonDividedDifference(u, sinu)
        # my_cos = NewtonDividedDifference(u, cosu)

        ################################################################################

        # h = CubicSpline(x, y, 'clamp')
        # ic(h(0.5))
        # h = CubicSpline(x, y, 'not-a-knot')
        # ic(h(0.5))
        # h = CubicSpline(x, y)
        # ic(h(0.5))
        # plt.figure(figsize=(8, 6))
        # t = np.linspace(0, 9, 500)
        # plt.plot(t, CubicSpline(x, y, 'natural')(t), label='natural')
        # plt.plot(t, CubicSpline(x, y, 'clamp')(t), label='clamp')
        # plt.plot(t, CubicSpline(x, y, 'not-a-knot')(t), label='not-a-knot')
        # plt.scatter(x, y)
        # plt.legend()
        # plt.show()

        # s = Bezier(x, y)
        # ic(s(0.5))
        # plt.figure(figsize=(8, 6))
        # t = np.linspace(0, 1, 500)
        # xs, ys = s(t)
        # plt.plot(x, y)
        # plt.plot(xs, ys)
        # plt.scatter(x, y)
        # plt.show()
        
        # plt.figure(figsize=(8, 6))
        # t = np.linspace(0, 9, 500)
        # plt.plot(t, f(t), label="Lagrange")
        # plt.plot(t, g(t), label="Newton")
        # plt.plot(t, h(t), label="Cubic")
        # plt.plot(*s(np.linspace(0, 1, 500)), label="Bezier")
        # plt.scatter(x, y)
        # plt.legend()
        # plt.show()