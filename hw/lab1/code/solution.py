import poisson
import torch
import torch.nn as nn

class Test(object):
    def __init__(self):
        #* constants
        self.N = 10
        self.N_1 = self.N - 1
        self.h = 1 / self.N
        self.epsilon = 1e-6

        #* precise solution
        x, y = torch.meshgrid(
            torch.linspace(0.0, 1.0, self.N + 1, dtype=torch.float64),
            torch.linspace(0.0, 1.0, self.N + 1, dtype=torch.float64),
            indexing="xy"
        )
        x = x[1:-1, 1:-1].cuda()
        y = y[1:-1, 1:-1].cuda()
        self.u_star = torch.sin(torch.pi * x) * torch.sin(torch.pi * y)
        self.u_star = nn.ZeroPad2d(1)(self.u_star)

        #* A, b
        self.A = poisson.A_mat(self.N_1)
        self.b = poisson.fun(x, y).T.reshape((self.N_1 * self.N_1, 1)) * self.h ** 2

class Test1(Test):
    def __problem(self, name: str, k: int, uh: torch.DoubleTensor, x1: torch.DoubleTensor):
        print("\033\13331mUsing {} iteration...\033\1330m".format(name))
        print("\033\13333mParameters: h = {}, N = {}, epsilon = 1e-6\033\1330m".format(self.h, self.N))
        print("\033\13335m||x1 - x0||: {}\033\1330m".format(torch.max(abs(x1))))
        print("\033\13332mIteration times: {}\033\1330m".format(k))
        print("\033\13334mError: {}\033\1330m\n".format(torch.max(torch.abs(uh - self.u_star))))

    def problem1(self):
        k, uh, x1 = poisson.Jacobi(self.A, self.b, self.N_1, self.epsilon)
        self.__problem("Jacobi", k, uh, x1)

    def problem2(self):
        k, uh, x1 = poisson.GaussSeidel(self.A, self.b, self.N_1, self.epsilon)
        self.__problem("GaussSeidel", k, uh, x1)

    def problem3(self):
        k, uh, x1 = poisson.SOR(self.A, self.b, 1.2, self.N_1, self.epsilon)
        self.__problem("SOR, omega = 1.2", k, uh, x1)

    def problem4(self):
        k, uh, x1 = poisson.SOR(self.A, self.b, 1.3, self.N_1, self.epsilon)
        self.__problem("SOR, omega = 1.3", k, uh, x1)

    def problem5(self):
        k, uh, x1 = poisson.SOR(self.A, self.b, 1.9, self.N_1, self.epsilon)
        self.__problem("SOR, omega = 1.9", k, uh, x1)

    def problem6(self):
        k, uh, x1 = poisson.SOR(self.A, self.b, 0.9, self.N_1, self.epsilon)
        self.__problem("SOR, omega = 0.9", k, uh, x1)

    def problem7(self):
        k, uh, x1 = poisson.SOR(self.A, self.b, 1.527864045, self.N_1, self.epsilon)
        self.__problem("SOR, omega = 1.527864045", k, uh, x1)


class Test2(Test):
    def change(self, N: int):
        #* constants
        self.N = N
        self.N_1 = self.N - 1
        self.h = 1 / self.N
        # self.epsilon = 1e-6

        #* precise solution
        x, y = torch.meshgrid(
            torch.linspace(0.0, 1.0, self.N + 1, dtype=torch.float64),
            torch.linspace(0.0, 1.0, self.N + 1, dtype=torch.float64),
            indexing="xy"
        )
        x = x[1:-1, 1:-1].cuda()
        y = y[1:-1, 1:-1].cuda()
        self.u_star = torch.sin(torch.pi * x) * torch.sin(torch.pi * y)
        self.u_star = nn.ZeroPad2d(1)(self.u_star)

        #* A, b
        self.A = poisson.A_mat(self.N_1)
        self.b = poisson.fun(x, y).T.reshape((self.N_1 * self.N_1, 1)) * self.h ** 2

    def __problem(self, N: int):
        self.change(N)

        print("\033\13331mUsing Jacobi iteration...\033\1330m")
        print("\033\13333mParameters: h = {}, N = {}, epsilon = 1e-6\033\1330m".format(self.h, self.N))

        k, uh, x1 = poisson.Jacobi(self.A, self.b, self.N_1, self.epsilon)

        print("\033\13335m||x1 - x0||: {}\033\1330m".format(torch.max(abs(x1))))
        print("\033\13332mIteration times: {}\033\1330m".format(k))
        print("\033\13334mError: {}\033\1330m\n".format(torch.max(torch.abs(uh - self.u_star))))

    def problem1(self):
        self.__problem(10)

    def problem2(self):
        self.__problem(20)

    def problem3(self):
        self.__problem(50)

    def problem4(self):
        self.__problem(100)

if __name__ == "__main__":

    test1 = Test1()
    test1.problem1()
    test1.problem2()
    test1.problem3()
    test1.problem4()
    test1.problem5()
    test1.problem6()
    test1.problem7()

    test2 = Test2()
    test2.problem1()
    test2.problem2()
    test2.problem3()
    test2.problem4()
