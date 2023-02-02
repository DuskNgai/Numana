from typing import Tuple

import numpy as np
import torch
import torch.nn as nn

def fun(x: torch.DoubleTensor, y: torch.DoubleTensor) -> torch.DoubleTensor:
    return 2 * torch.pi * torch.pi * torch.sin(torch.pi * x) * torch.sin(torch.pi * y)

def A_mat(N_1: int) -> torch.DoubleTensor:
    # C is a matrix where element on its -1 and 1 diagonal are 1s,
    # others are 0s
    C = np.eye(N_1, k = 1, dtype=np.float64) + np.eye(N_1, k = -1, dtype=np.float64)
    C = torch.DoubleTensor(C).cuda()
    # M is a matrix where element on its main diagonal are 4s,
    # -1 and 1 diagonal are 1s, others are 0s
    M = 4 * torch.eye(N_1, dtype=torch.float64).cuda() - C

    A_size = N_1 * N_1
    A = -np.eye(A_size, k = N_1, dtype=np.float64) - np.eye(A_size, k = -N_1, dtype=np.float64)
    A = torch.DoubleTensor(A).cuda()
    for i in range(N_1):
        A[i * N_1: (i + 1) * N_1, i * N_1: (i + 1) * N_1] = M

    return A

def iteration(
    B: torch.DoubleTensor, f: torch.DoubleTensor, N_1: int, epsilon: float
) -> Tuple[int, torch.DoubleTensor]:
    print("\033\13331mSpectrum radius = {}\033\1330m".format(torch.max(torch.abs(torch.linalg.eigvals(B)))))

    A_size = N_1 * N_1
    k = 0
    uh = torch.zeros((A_size, 1), dtype=torch.float64).cuda()
    while True:
        uh_prev = uh
        uh = B @ uh_prev + f
        norm_infty = torch.max(torch.abs(uh - uh_prev))
        if norm_infty < epsilon:
            break
        k += 1

    uh = nn.ZeroPad2d(1)(uh.reshape((N_1, N_1)))
    return k, uh, f

def Jacobi(
    A: torch.DoubleTensor, b: torch.DoubleTensor, N_1: int, epsilon: float
) -> Tuple[int, torch.DoubleTensor]:
    L = -torch.tril(A, diagonal=-1)
    U = -torch.triu(A, diagonal= 1)
    J = 0.25 * (L + U)
    f = 0.25 * b

    return iteration(J, f, N_1, epsilon)

def GaussSeidel(
    A: torch.DoubleTensor, b: torch.DoubleTensor, N_1: int, epsilon: float
) -> Tuple[int, torch.DoubleTensor]:
    A_size = N_1 * N_1
    L = -torch.tril(A, diagonal=-1)
    U = -torch.triu(A, diagonal= 1)
    D = 4.0 * torch.eye(A_size, dtype=torch.float64).cuda()
    DL_inv = torch.linalg.inv(D - L)
    G = DL_inv @ U
    f = DL_inv @ b

    return iteration(G, f, N_1, epsilon)

def SOR(
    A: torch.DoubleTensor, b: torch.DoubleTensor, omega: float, N_1: int, epsilon: float
) -> Tuple[int, torch.DoubleTensor]:
    A_size = N_1 * N_1
    L = -torch.tril(A, diagonal=-1)
    U = -torch.triu(A, diagonal= 1)
    D = 4.0 * torch.eye(A_size, dtype=torch.float64).cuda()
    D_omegaL_inv = torch.linalg.inv(D - omega * L)
    S = D_omegaL_inv @ ((1 - omega) * D + omega * U)
    f = omega * D_omegaL_inv @ b

    return iteration(S, f, N_1, epsilon)
