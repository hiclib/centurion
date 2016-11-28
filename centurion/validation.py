import numpy as np


def compute_error(results, groundtruth):
    """
    Computes error

    Parameters
    ----------
    results : ndarray (L, )

    groundtruth : ndarray (L, 2)

    Returns
    -------
    error : ndarray (L, )
    """
    mask = np.isnan(groundtruth[:, 0])
    error = (results[:, np.newaxis] - groundtruth)
    error[error[:, 0] > 0, 0] = 0
    error[error[:, 1] < 0, 1] = 0
    error = error.sum(axis=1)
    if np.any(mask):
        error[mask] = np.nan
    return error
