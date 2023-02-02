import numpy as np
import matplotlib.pyplot as plt

import ode

plt.rcParams['font.family'] = ["Georgia"]
plt.rcParams['font.size'] = 24.0

class Test(object):
    def draw(self, xs: np.ndarray, ys: np.ndarray, title: str, filename: str):
        plt.figure(figsize=(10, 8))
        plt.plot(xs, ys, label="Numerical", color="orange")
        plt.plot(self.x_real, self.y_real, label="Symbolic", color="purple")
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig(filename)
        plt.close()

class Test1(Test):
    def __init__(self):
        self.initial = 100.0
        self.interval = (0.0, 1.0)
        self.x_real = np.linspace(*self.interval, 100)
        self.y_real = self.initial * np.exp(-50 * np.linspace(*self.interval, 100))

    def problem1(self):
        print("\033\13331mUsing Explicit Euler's Methods...\033\1330m")
        self.solver = ode.ExplicitEuler(lambda x, y: -50 * y)
        for h in range(5, 0, -1):
            step_size = h * 0.01
            xs, ys = self.solver.solve(self.initial, step_size, self.interval)
            self.draw(
                xs, ys,
                r"Explicit Euler's Methods with $h={:.2f}$".format(step_size),
                r"lab2\figure\Q11_h_{:.2f}.png".format(step_size)
            )

    def problem2(self):
        print("\033\13331mUsing Direct Implicit Euler's Methods...\033\1330m")
        self.solver = ode.ImplicitEuler(lambda h, x, y: y / (1 + 50 * h))
        for h in range(5, 0, -1):
            step_size = h * 0.01
            xs, ys = self.solver.solve(self.initial, step_size, self.interval)
            self.draw(
                xs, ys,
                r"Implicit Euler's Methods with $h={:.2f}$".format(step_size),
                r"lab2\figure\Q12_h_{:.2f}_precise.png".format(step_size)
            )

    def problem3(self):
        print("\033\13331mUsing Iterative Implicit Euler's Methods...\033\1330m")
        self.solver = ode.ImplicitEuler(lambda x, y: -50 * y)
        for h in range(5, 0, -1):
            step_size = h * 0.01
            xs, ys = self.solver.iterativeSolve(self.initial, step_size, self.interval)
            self.draw(
                xs, ys,
                r"Implicit Euler's Methods with $h={:.2f}$".format(step_size),
                r"lab2\figure\Q12_h_{:.2f}.png".format(step_size)
            )

class Test2(Test):
    def __init__(self):
        self.initial = [3.0, 2.0]
        self.interval = (0.0, 4 * np.pi)
        self.x_real = np.cos(np.linspace(0.0, 2.0 * np.pi, 100)) * (3 * np.sqrt(2))
        self.y_real = np.sin(np.linspace(0.0, 2.0 * np.pi, 100)) * (2 * np.sqrt(2))

    def problem1(self):
        print("\033\13331mUsing Trapezoid Euler's Methods...\033\1330m")
        self.solver = ode.TrapezoidEulerSystem(
            lambda x, y: 4.5 * y,
            lambda x, y: -2.0 * x,
        )
        for exp in range(-1, -4, -1):
            h = 10 ** exp
            xs, ys = self.solver.iterativeSolve(self.initial, h, self.interval)
            self.draw(
                xs, ys,
                r"Trapezoid Euler's Method with $h={:.3f}$".format(h),
                r"lab2\figure\Q21_trapezoid_h_{:.3f}.png".format(h)
            )

    def problem2(self):
        print("\033\13331mUsing Explicit Euler's Methods...\033\1330m")
        self.solver = ode.ExplicitEulerSystem(
            lambda x, y: 4.5 * y,
            lambda x, y: -2.0 * x,
        )
        for exp in range(-1, -4, -1):
            h = 10 ** exp
            xs, ys = self.solver.solve(self.initial, h, self.interval)
            self.draw(
                xs, ys,
                r"Explicit Euler's Method with $h={:.3f}$".format(h),
                r"lab2\figure\Q22_explicit_h_{:.3f}.png".format(h)
            )

    def problem3(self):
        print("\033\13331mUsing Implicit Euler's Methods...\033\1330m")
        self.solver = ode.ImplicitEulerSystem(
            lambda x, y: 4.5 * y,
            lambda x, y: -2.0 * x,
        )
        for exp in range(-1, -4, -1):
            h = 10 ** exp
            xs, ys = self.solver.iterativeSolve(self.initial, h, self.interval)
            self.draw(
                xs, ys,
                r"Implicit Euler's Method with $h={:.3f}$".format(h),
                r"lab2\figure\Q23_implicit_h_{:.3f}.png".format(h)
            )

    def problem4(self):
        print("\033\13331mUsing Forth Order of Runge Kutta...\033\1330m")
        self.solver = ode.RungeKuttaSystem(
            lambda x, y: 4.5 * y,
            lambda x, y: -2.0 * x,
        )
        for exp in range(-1, -4, -1):
            h = 10 ** exp
            xs, ys = self.solver.solve(self.initial, h, self.interval)
            self.draw(
                xs, ys,
                r"Runge Kutta with $h={:.3f}$".format(h),
                r"lab2\figure\Q24_rk4_h_{:.3f}.png".format(h)
            )

class Test3(Test):
    def __init__(self):
        self.initial = [np.pi / 3, -0.5]
        self.interval = [0.0, 10.0]

    def problem1(self):
        print("\033\13331mUsing Implicit Euler's Methods...\033\1330m")
        print("\033\13331mUsing Trapezoid Euler's Methods...\033\1330m")
        h = 0.02
        solver1 = ode.ImplicitEulerSystem(
            lambda y1, y2: y2,
            lambda y1, y2: -np.sin(y1),
        )
        solver2 = ode.TrapezoidEulerSystem(
            lambda y1, y2: y2,
            lambda y1, y2: -np.sin(y1),
        )

        yi1, yi2 = solver1.iterativeSolve(self.initial, h, self.interval)
        yt1, yt2 = solver2.iterativeSolve(self.initial, h, self.interval)
        xs = np.linspace(*self.interval, yi1.shape[0])
        plt.figure(figsize=(10, 8))
        plt.plot(xs, yi1, label="Implicit", color="orange")
        plt.plot(xs, yt1, label="Trapezoid", color="purple")
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title(r"Different Methods with $h={:.2f}$".format(h))
        plt.legend()
        plt.tight_layout()
        plt.savefig(r"lab2\figure\Q31_h_{:.2f}.png".format(h))
        plt.close()

    def problem2(self):
        print("\033\13331mUsing Implicit Euler's Methods...\033\1330m")
        print("\033\13331mUsing Trapezoid Euler's Methods...\033\1330m")
        h = 0.02
        solver1 = ode.ImplicitEulerSystem(
            lambda y1, y2: y2,
            lambda y1, y2: -np.sin(y1),
        )
        solver2 = ode.TrapezoidEulerSystem(
            lambda y1, y2: y2,
            lambda y1, y2: -np.sin(y1),
        )
        solver3 = ode.RungeKuttaSystem(
            lambda y1, y2: y2,
            lambda y1, y2: -np.sin(y1),
        )
        yi1, yi2 = solver1.iterativeSolve(self.initial, h, self.interval)
        yt1, yt2 = solver2.iterativeSolve(self.initial, h, self.interval)
        yr1, yr2 = solver3.solve(self.initial, h, self.interval)
        xs = np.linspace(*self.interval, yi1.shape[0])
        plt.figure(figsize=(10, 8))
        plt.plot(xs, yi1, label="Implicit", color="orange")
        plt.plot(xs, yt1, label="Trapezoid", color="purple")
        plt.plot(xs, yr1, label="RK4", color="pink")
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title(r"Different Methods with $h={:.2f}$".format(h))
        plt.legend()
        plt.tight_layout()
        plt.savefig(r"lab2\figure\Q32_h_{:.2f}.png".format(h))
        plt.close()

if __name__ == "__main__":

    test1 = Test1()
    test1.problem1()
    test1.problem2()
    test1.problem3()

    test2 = Test2()
    test2.problem1()
    test2.problem2()
    test2.problem3()
    test2.problem4()

    test3 = Test3()
    test3.problem1()
    test3.problem2()
