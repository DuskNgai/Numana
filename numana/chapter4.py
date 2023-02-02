import numpy as np

class LeastSquare(object):

    @staticmethod
    def solve(A: np.ndarray, b: np.ndarray):
        """
        Find the least square solution of the equation `Ax = b`.
        @param `A, b`: The array to solve.
        @return: The least square solution.
        @return: The RMS error of the least square solution.
        """

        assert len(A.shape) == 2, "A must be a 2D array."
        # A.size() = (M, N)
        # x.size() = (N, P)
        # b.size() = (M, P)

        AT = A.T
        ATA = AT @ A
        ATb = AT @ b

        x = np.linalg.inv(ATA) @ ATb
        err = np.sqrt(np.mean(np.square(b - A @ x), axis=0))
        return x, err
