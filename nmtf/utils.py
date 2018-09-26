import numpy as np


def dispersion_coefficient(X):
    coeff = 0
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            coeff += 4*(X[i,j] - 0.5)**2
    return coeff/X.shape[0]**2


def connectivity_matrix(X):
    indices = np.argmax(X, axis=1)
    C = np.zeros_like(X)
    for r, i in enumerate(indices):
        C[r,i] = 1
    return C.dot(C.T)


def squared_norm(x):
    """Squared Euclidean or Frobenius norm of x.
    Faster than norm(x) ** 2.
    Parameters
    ----------
    x : array_like
    Returns
    -------
    float
        The Euclidean norm when x is a vector, the Frobenius norm when x
        is a matrix (2-d array).
    """
    x = np.ravel(x, order='K')
    if np.issubdtype(x.dtype, np.integer):
        warnings.warn('Array type is integer, np.dot may overflow. '
                      'Data should be float type to avoid this issue',
                      UserWarning)
    return np.dot(x, x)
