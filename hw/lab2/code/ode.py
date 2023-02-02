from typing import Callable, Sequence, Tuple

import numpy as np

#! For single ode
class Solver(object):
    def __init__(self, f: Callable[[Sequence[float]], float]) -> None:
        self.f = f
        self.num_iter = 9

    def check(self, initial: float, step_size: float, interval: Tuple[float, float]) -> None:
        """Check validity of the variables."""
        start, end = interval
        assert isinstance(initial, (int, float)), "Initial must be a number."
        assert isinstance(step_size, float) and 0 < step_size < 1, "Step_size should be a small number."
        assert isinstance(start, (int, float)) and isinstance(end, (int, float)) and start < end, "Invalid interval."

class ExplicitEuler(Solver):
    def solve(self, initial: float, step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> Tuple[np.ndarray, ...]:
        self.check(initial, step_size, interval)
        start, end = interval

        #* y_{n + 1} = y_{n} + h * f(x_{n}, y_{n})
        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = np.zeros_like(xs)
        ys[0] = initial
        for i in range(xs.shape[0] - 1):
            ys[i + 1] = ys[i] + self.f(xs[i], ys[i]) * step_size

        return xs, ys

class ImplicitEuler(Solver):
    def solve(self, initial: float, step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> Tuple[np.ndarray, ...]:
        self.check(initial, step_size, interval)
        start, end = interval

        #* y_{n + 1} = y_{n} + h * f(x_{n + 1}, y_{n + 1})
        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = np.zeros_like(xs)
        ys[0] = initial
        for i in range(xs.shape[0] - 1):
            ys[i + 1] = self.f(step_size, xs[i], ys[i])

        return xs, ys

    def iterativeSolve(self, initial: float, step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> Tuple[np.ndarray, ...]:
        self.check(initial, step_size, interval)
        start, end = interval

        #* y_{n + 1} = y_{n} + h * f(x_{n + 1}, y_{n + 1})
        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = np.zeros_like(xs)
        ys[0] = initial
        for i in range(xs.shape[0] - 1):
            #* using explicit Euler's method for initial value
            y_i_prev = ys[i] + self.f(xs[i], ys[i]) * step_size
            y_i_next = ys[i] + self.f(xs[i + 1], y_i_prev) * step_size

            for _ in range(self.num_iter):
                y_i_prev = y_i_next
                y_i_next = ys[i] + self.f(xs[i + 1], y_i_prev) * step_size
                if np.allclose(y_i_prev, y_i_next):
                    break
            ys[i + 1] = y_i_next

        return xs, ys

#! For ode system
class SystemSolver(object):
    def __init__(self, *f: Sequence[Callable[[Sequence[float]], float]]) -> None:
        self.f = f
        self.num_iter = 99

    def check(self, initial: list[float], step_size: float, interval: Tuple[float, float]) -> None:
        """Check validity of the variables."""
        start, end = interval
        for init in initial:
            assert isinstance(init, (int, float)), "Initial must be a number."
        assert isinstance(step_size, float) and 0 < step_size < 1, "Step_size should be a small number."
        assert isinstance(start, (int, float)) and isinstance(end, (int, float)) and start < end, "Invalid interval."

class ExplicitEulerSystem(SystemSolver):
    def solve(self, initial: list[float], step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> list[np.ndarray]:
        self.check(initial, step_size, interval)
        start, end = interval

        #* y_{1, n + 1} = y_{1, n} + h * f_{1}(y_{1, n}, y_{2, n}, ...)
        #* y_{2, n + 1} = y_{2, n} + h * f_{2}(y_{1, n}, y_{2, n}, ...)
        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = []
        for init in initial:
            y = np.zeros_like(xs)
            y[0] = init
            ys.append(y)

        #* y_{i, j + 1} = y_{i, j} + h * f_{i}(y_{1, j}, y_{2, j}, ...)
        for j in range(xs.shape[0] - 1):
            y_i_j = [y_i[j] for y_i in ys]
            for y_i, f_i in zip(ys, self.f):
                y_i[j + 1] = y_i[j] + f_i(*y_i_j) * step_size

        return ys

class ImplicitEulerSystem(SystemSolver):
    def iterativeSolve(self, initial: list[float], step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> list[np.ndarray]:
        self.check(initial, step_size, interval)
        start, end = interval

        #* y_{1, n + 1} = y_{1, n} + h * f_{i}(y_{1, n + 1}, y_{2, n + 1}, ...)
        #* y_{2, n + 1} = y_{2, n} + h * f_{i}(y_{1, n + 1}, y_{2, n + 1}, ...)
        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = []
        for init in initial:
            y = np.zeros_like(xs)
            y[0] = init
            ys.append(y)

        #* y_{i, j + 1} = y_{i, j} + h * f_{i}(y_{1, j + 1}, y_{2, j + 1}, ...)
        for j in range(xs.shape[0] - 1):
            #* using explicit Euler's method for initial value
            y_i_j = [y_i[j] for y_i in ys]
            y_i_j_prev = [y_i[j] + f_i(*y_i_j) * step_size for y_i, f_i in zip(ys, self.f)]
            y_i_j_next = [y_i[j] + f_i(*y_i_j_prev) * step_size for y_i, f_i in zip(ys, self.f)]

            for _ in range(self.num_iter):
                y_i_j_prev = y_i_j_next
                y_i_j_next = [y_i[j] + f_i(*y_i_j_prev) * step_size for y_i, f_i in zip(ys, self.f)]
                if np.allclose(y_i_j_prev, y_i_j_next):
                    break

            for y_i, y_i_j1 in zip(ys, y_i_j_next):
                y_i[j + 1] = y_i_j1

        return ys

class TrapezoidEulerSystem(SystemSolver):
    def iterativeSolve(self, initial: list[float], step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> list[np.ndarray]:
        self.check(initial, step_size, interval)
        start, end = interval

        #* y_{1, n + 1} = y_{1, n} + h * [f_{1}(y_{1, n}, y_{2, n}, ...) + f_{i}(y_{1, n + 1}, y_{2, n + 1}, ...)]
        #* y_{2, n + 1} = y_{2, n} + h * [f_{2}(y_{1, n}, y_{2, n}, ...) + f_{i}(y_{1, n + 1}, y_{2, n + 1}, ...)]
        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = []
        for init in initial:
            y = np.zeros_like(xs)
            y[0] = init
            ys.append(y)

        #* y_{i, j + 1} = y_{i, j} + h * 0.5 * [f_{i}(y_{1, j}, y_{2, j}, ...) + f_{i}(y_{1, j + 1}, y_{2, j + 1}, ...)]
        for j in range(xs.shape[0] - 1):
            #* using explicit Euler's method for initial value
            y_i_j = [y_i[j] for y_i in ys]
            y_i_j_prev = [y_i[j] + f_i(*y_i_j) * step_size for y_i, f_i in zip(ys, self.f)]
            y_i_j_next = [y_i[j] + f_i(*y_i_j_prev) * step_size for y_i, f_i in zip(ys, self.f)]

            for _ in range(self.num_iter):
                y_i_j_prev = y_i_j_next
                y_i_j_next = [y_i[j] + f_i(*y_i_j_prev) * step_size for y_i, f_i in zip(ys, self.f)]
                if np.allclose(y_i_j_prev, y_i_j_next):
                    break

            y_i_j = [y_i[j] for y_i in ys]
            f_y_i_j = [f_i(*y_i_j) for f_i in self.f]
            f_y_i_j1 = [f_i(*y_i_j_next) for f_i in self.f]
            for y_i, a, b in zip(ys, f_y_i_j, f_y_i_j1):
                y_i[j + 1] = y_i[j] + (a + b) * step_size * 0.5

        return ys

class RungeKuttaSystem(SystemSolver):
    def rk4(self, step_size: float, y_i_j: Sequence[float]) -> float:
        K1 = [f_i(*y_i_j) for f_i in self.f]
        y1 = [y_j + k1_j * step_size * 0.5 for y_j, k1_j in zip(y_i_j, K1)]

        K2 = [f_i(*y1) for f_i in self.f]
        y2 = [y_j + k2_j * step_size * 0.5 for y_j, k2_j in zip(y_i_j, K2)]

        K3 = [f_i(*y2) for f_i in self.f]
        y3 = [y_j + k2_j * step_size for y_j, k2_j in zip(y_i_j, K3)]

        K4 = [f_i(*y3) for f_i in self.f]

        y_i_j_next = [y_j + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * step_size / 6.0 for y_j, k1, k2, k3, k4 in zip(y_i_j, K1, K2, K3, K4)]
        return y_i_j_next

    def solve(self, initial: list[float], step_size: float, interval: Tuple[float, float] = (0.0, 1.0)) -> list[np.ndarray]:
        self.check(initial, step_size, interval)
        start, end = interval

        xs = np.arange(start, end + step_size, step_size)
        xs[-1] = end
        ys = []
        for init in initial:
            y = np.zeros_like(xs)
            y[0] = init
            ys.append(y)

        for j in range(xs.shape[0] - 1):
            y_i_j = [y_i[j] for y_i in ys]
            y_i_j_next = self.rk4(step_size, y_i_j)

            for y_i, y_i_j1 in zip(ys, y_i_j_next):
                y_i[j + 1] = y_i_j1

        return ys
