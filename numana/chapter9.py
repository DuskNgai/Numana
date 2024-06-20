import math
import random
from typing import Callable

import numpy as np

class RandomWalk(object):
    @staticmethod
    def reach_top_endpoint(endpoints: tuple[int, int], prob: float) -> int:
        x = 0
        num_success = 0
        while True:
            if random.random() < prob:
                x += 1
            else:
                x -= 1

            if x == endpoints[0]:
                break
            elif x == endpoints[1]:
                num_success += 1
                break

        return num_success

    @staticmethod
    def num_steps_to_endpoint(endpoints: tuple[int, int], prob: float) -> int:
        x = 0
        num_steps = 0
        while x not in endpoints:
            if random.random() < prob:
                x += 1
            else:
                x -= 1

            num_steps += 1

        return num_steps

class Solver(object):
    """
    Solving the Stochastic Differential Equation (SDE) of the form:
    `dX(t) = f(X(t), t) dt + g(X(t), t) dw(t)` with the numerical method.
    """

    def __init__(self,
        f: Callable[[float, float, float], float],
        g: Callable[[float, float, float], float],
        x_hat: Callable[[float, float, float], float],
    ) -> None:
        self.f = f
        self.g = g
        self.x_hat = x_hat

    @staticmethod
    def wiener_process(t0: float, t1: float, n: int) -> np.ndarray:
        """
        Generate the Wiener process `w(t)` with the number of steps `n` in the interval `[t0, t1]`.
        `w(0) = 0` and `w(t_{i+1}) - w(t_i) ~ N(0, t_{i+1} - t_i)`.

        Returns:
            (np.ndarray): The Wiener process `w(t)` with shape [n + 1].
        """
        ts = np.linspace(t0, t1, n + 1, dtype=float)
        dts = np.diff(ts)
        dws = np.sqrt(dts) * np.random.randn(n)
        ws = np.insert(np.cumsum(dws), 0, 0)
        return ws
    
    def accurate_solution(self, x0: float, t0: float, t1: float, n: int, ws: np.ndarray | None) -> np.ndarray:
        """
        Compute the accurate solution of the SDE.

        Returns:
            (np.ndarray): The accurate solution of the SDE with shape [n + 1].
        """
        if ws is None:
            ws = self.wiener_process(t0, t1, n)

        ts = np.linspace(t0, t1, n + 1, dtype=float)
        xs = self.x_hat(ws, x0, ts)
        return np.array(xs)

    def euler_maruyama(self, x0: float, t0: float, t1: float, n: int, ws: np.ndarray | None) -> np.ndarray:
        """
        Compute the solution of the SDE by the Euler-Maruyama method.
        """
        if ws is None:
            ws = self.wiener_process(t0, t1, n) # [n + 1]

        ts = np.linspace(t0, t1, n + 1, dtype=float)
        dts = np.diff(ts)
        dws = np.diff(ws)

        x, xs = x0, [x0]
        for w, dw, t, dt in zip(ws[:-1], dws, ts[:-1], dts):
            x += self.f(w, x, t) * dt + self.g(w, x, t) * dw
            xs.append(x)

        return np.array(xs)

    def milstein(self, g_x: Callable[[float, float], float], x0: float, t0: float, t1: float, n: int, ws: np.ndarray | None) -> np.ndarray:
        """
        Compute the solution of the SDE by the Milstein method.
        """
        if ws is None:
            ws = self.wiener_process(t0, t1, n)

        ts = np.linspace(t0, t1, n + 1, dtype=float)
        dts = np.diff(ts)
        dws = np.diff(ws)

        x, xs = x0, [x0]
        for w, dw, t, dt in zip(ws[:-1], dws, ts[:-1], dts):
            x += self.f(w, x, t) * dt + self.g(w, x, t) * dw + 0.5 * self.g(w, x, t) * g_x(w, x, t) * (dw ** 2 - dt)
            xs.append(x)

        return np.array(xs)

    def runge_kutta(self, x0: float, t0: float, t1: float, n: int, ws: np.ndarray | None) -> np.ndarray:
        """
        Compute the solution of the SDE by the Runge-Kutta method.
        """
        if ws is None:
            ws = self.wiener_process(t0, t1, n)

        ts = np.linspace(t0, t1, n + 1, dtype=float)
        dts = np.diff(ts)
        dws = np.diff(ws)

        x, xs = x0, [x0]
        for w, dw, t, dt in zip(ws[:-1], dws, ts[:-1], dts):
            dt_sqrt = np.sqrt(dt)
            g_x_t = self.g(w, x, t)
            x += self.f(w, x, t) * dt + g_x_t * dw + 0.5 * (self.g(w, x + g_x_t * dt_sqrt, t) - g_x_t) / dt_sqrt * (dw ** 2 - dt)
            xs.append(x)

        return np.array(xs)
