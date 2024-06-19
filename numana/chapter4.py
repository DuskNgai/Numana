import numpy as np

class LeastSquare(object):

    @staticmethod
    def solve(A: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Find the least square solution of the equation `Ax = b`.

        Args:
            A (np.ndarray): The matrix A.
            b (np.ndarray): The right hand side of the equation.

        Returns:
            Tuple[np.ndarray, np.ndarray]: The solution x and the error.
        """

        assert A.ndim == 2, "A must be a 2D array."
        # A.size() = (M, N)
        # x.size() = (N, P)
        # b.size() = (M, P)

        AT = A.T
        ATA = AT @ A
        ATb = AT @ b

        x = np.linalg.inv(ATA) @ ATb
        err = np.sqrt(np.mean(np.square(b - A @ x), axis=0))
        return x, err
