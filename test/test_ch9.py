import random
from pathlib import Path
import sys
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pytest
from rich import print

sys.path.append(str(Path(__file__).resolve().parents[1]))

import numana.chapter9 as ch9

class TestRandomWalk(object):
    @staticmethod
    def seed(seed: int) -> None:
        random.seed(seed)
        np.random.seed(seed)

    def output_reach_top_endpoint(self, endpoints: tuple[int, int], prob: float) -> int:
        TestRandomWalk.seed(0)

        num_success = 0
        for _ in range(self.num_trials):
            num_success += ch9.RandomWalk.reach_top_endpoint(endpoints, prob)

        b, a = endpoints
        if prob == 0.5:
            expected_value = - b / (a - b)
        else:
            p, q = prob, 1 - prob
            expected_value = ((q / p) ** (-b) - 1) / ((q / p) ** (a - b) - 1)
        print("The frequency of reaching the top endpoint of {} with probability {} is {}, while the expected value is {}.".format(endpoints, prob, num_success / self.num_trials, expected_value))

    def test_reach_top_endpoint(self) -> None:
        self.num_trials = 10000

        # (1)
        self.output_reach_top_endpoint((-2, 5), 0.5)
        self.output_reach_top_endpoint((-5, 3), 0.5)
        self.output_reach_top_endpoint((-8, 3), 0.5)

        # (3)
        self.output_reach_top_endpoint((-2, 5), 0.7)
        self.output_reach_top_endpoint((-5, 3), 0.7)
        self.output_reach_top_endpoint((-8, 3), 0.7)

    def output_num_steps_to_endpoint(self, endpoints: tuple[int, int], prob: float) -> int:
        TestRandomWalk.seed(0)

        num_steps = 0
        for _ in range(self.num_trials):
            num_steps += ch9.RandomWalk.num_steps_to_endpoint(endpoints, prob)

        b, a = endpoints
        if prob == 0.5:
            expected_value = - a * b
        else:
            p, q = prob, 1 - prob
            expected_value = (-b - (a - b) * (1 - (q / p) ** (-b)) / (1 - (q / p) ** (a - b))) / (q - p)
        print("The average number of steps to reach the endpoint of {} with probability {} is {}, while the expected value is {}.".format(endpoints, prob, num_steps / self.num_trials, expected_value))

    def test_num_steps_to_endpoint(self) -> None:
        self.num_trials = 10000

        # (2)
        self.output_num_steps_to_endpoint((-2, 5), 0.5)
        self.output_num_steps_to_endpoint((-5, 3), 0.5)
        self.output_num_steps_to_endpoint((-8, 3), 0.5)

        # (4)
        self.output_num_steps_to_endpoint((-2, 5), 0.7)
        self.output_num_steps_to_endpoint((-5, 3), 0.7)
        self.output_num_steps_to_endpoint((-8, 3), 0.7)

class TestSolver(object):
    @staticmethod
    def seed(seed: int) -> None:
        random.seed(seed)
        np.random.seed(seed)

    def output_solver(self,
        f: Callable[[float, float, float], float],
        g: Callable[[float, float, float], float],
        g_x: Callable[[float, float, float], float],
        x_hat: Callable[[float, float, float], float],
        x0: float,
        t0: float,
        t1: float,
        dt: float
     ) -> float:
        TestSolver.seed(0)
        self.num_trials = 100

        solver = ch9.Solver(f, g, x_hat)
        n = int((t1 - t0) / dt)

        ws = solver.wiener_process(t0, t1, n)
        x1_euler_maruyama = np.array(solver.euler_maruyama(x0, t0, t1, n, ws))
        x1_milstein = np.array(solver.milstein(g_x, x0, t0, t1, n, ws))
        x1_runge_kutta = np.array(solver.runge_kutta(x0, t0, t1, n, ws))
        x1_hat = np.array(solver.accurate_solution(x0, t0, t1, n, ws))

        plt.figure(figsize=(8, 8))
        plt.plot(np.linspace(t0, t1, n + 1), x1_euler_maruyama, label="Euler-Maruyama")
        plt.plot(np.linspace(t0, t1, n + 1), x1_milstein, label="Milstein")
        plt.plot(np.linspace(t0, t1, n + 1), x1_runge_kutta, label="Runge-Kutta")
        plt.plot(np.linspace(t0, t1, n + 1), x1_hat, label="True")
        plt.grid()
        plt.legend()
        plt.show()

        print("The average error of the Euler-Maruyama method at time {} with the number of steps {} is {}.".format(t1, n, np.mean(np.abs(x1_euler_maruyama - x1_hat))))
        print("The average error of the Milstein method at time {} with the number of steps {} is {}.".format(t1, n, np.mean(np.abs(x1_milstein - x1_hat))))
        print("The average error of the Runge-Kutta method at time {} with the number of steps {} is {}.".format(t1, n, np.mean(np.abs(x1_runge_kutta - x1_hat))))
        print()

    def test_solver(self) -> None:
        # (0), Wiener process
        f = lambda w, x, t: 0
        g = lambda w, x, t: 1
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: w + x0
        x0, t0, t1, dt = 1, 0, 1, 0.1
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (0), Geometric Wiener process
        r, s = 0.1, 0.3
        f = lambda w, x, t: r * x
        g = lambda w, x, t: s * x
        g_x = lambda w, x, t: s
        x_hat = lambda w, x0, t: x0 * np.exp((r - 0.5 * s ** 2) * t + s * w)
        x0, t0, t1, dt = 1, 0, 2, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (0), Brownian bridge
        f = lambda w, x, t: (5 - x) / (5 - t)
        g = lambda w, x, t: 1
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: (5 - 1) / (5 - 1) * (t - 1) + 1
        x0, t0, t1, dt = 1, 1, 5, 0.001
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (1)
        f = lambda w, x, t: w
        g = lambda w, x, t: t
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: t * w + x0
        x0, t0, t1, dt = 0, 0, 10, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (2)
        f = lambda w, x, t: 0
        g = lambda w, x, t: 2 * w
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: w ** 2 - t + x0
        x0, t0, t1, dt = 0, 0, 10, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (3)
        f = lambda w, x, t: (1 - w ** 2) * np.exp(-2 * x)
        g = lambda w, x, t: 2 * w * np.exp(-x)
        g_x = lambda w, x, t: -2 * w * x * np.exp(-x)
        x_hat = lambda w, x0, t: np.log(1 + w ** 2) + x0
        x0, t0, t1, dt = 1, 0, 1, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (4)
        f = lambda w, x, t: w
        g = lambda w, x, t: np.cbrt(9 * x ** 2)
        g_x = lambda w, x, t: 2 / np.cbrt(3 * x)
        x_hat = lambda w, x0, t: w ** 3 / 3 + x0
        x0, t0, t1, dt = 1, 0, 1, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (5)
        f = lambda w, x, t: x * t
        g = lambda w, x, t: np.exp(t ** 2 / 2)
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: (1 + w) * np.exp(t ** 2 / 2)
        x0, t0, t1, dt = 1, 0, 2, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (6)
        f = lambda w, x, t: 0
        g = lambda w, x, t: 3 * (w ** 2 - t)
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: w ** 3 - 3 * t * w + x0
        x0, t0, t1, dt = 0, 0, 2, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (7), nan for runge_kutta
        f = lambda w, x, t: -0.5 * x
        g = lambda w, x, t: np.sqrt(1 - x ** 2)
        g_x = lambda w, x, t: -x / np.sqrt(1 - x ** 2)
        x_hat = lambda w, x0, t: np.sin(w) + x0
        x0, t0, t1, dt = 0, 0, 1, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (8)
        f = lambda w, x, t: x * (1 + 2 * np.log(x))
        g = lambda w, x, t: 2 * w * x
        g_x = lambda w, x, t: 2 * w
        x_hat = lambda w, x0, t: x0 * np.exp(w ** 2)
        x0, t0, t1, dt = 1, 0, 1, 0.01
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

        # (12)
        f = lambda w, x, t: -x
        g = lambda w, x, t: 1
        g_x = lambda w, x, t: 0
        x_hat = lambda w, x0, t: x0 * np.exp(-t) + w
        x0, t0, t1, dt = np.e, 0, 1, 0.001
        self.output_solver(f, g, g_x, x_hat, x0, t0, t1, dt)

if __name__ == "__main__":
    pytest.main(["-s", "test_ch9.py::TestRandomWalk::test_reach_top_endpoint"])
    pytest.main(["-s", "test_ch9.py::TestRandomWalk::test_num_steps_to_endpoint"])
    pytest.main(["-s", "test_ch9.py::TestSolver::test_solver"])
