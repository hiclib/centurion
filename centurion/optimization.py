import numpy as np
from scipy import optimize


def objective_function(x, lengths=None, counts=None):
    """
    Objective function for the centromere position
    """

    p = x[:len(lengths)]
    a = x[len(lengths)]
    b = x[len(lengths) + 1]
    sigma = x[len(lengths) + 2]

    p1 = np.array([i for i, l in enumerate(lengths) for j in range(int(l))])
    begin = np.concatenate([[0], lengths.cumsum()])[:-1]
    end = lengths.cumsum() - 1
    penalty = 0
    # FIXME this is not useful. Should add proper constraints
    if np.any(p[p > end]) or np.any(p[p < begin]):
        d = ((p[p > end] - end[p > end]) ** 2).sum()
        d += ((p[p < begin] - begin[p < begin]) ** 2).sum()
        penalty = 10e8 * d

    if counts is None:
        counts = 1

    return lambda x, y: (counts != 0) * (
        a * np.exp(
            - (p[p1[x]] - x) ** 2 / (2 * sigma ** 2)
            - (p[p1[y]] - y) ** 2 / (2 * sigma ** 2)) + b) + penalty


def fit_gaussian(data, init, lengths, counts=False, factor=1, ftol=1e-10,
                 gtol=1e-10, xtol=1e-10):
    if counts:
        counts = data
    else:
        counts = None
    errorfunction = lambda p: np.ravel(
        objective_function(p, lengths=lengths, counts=counts)(
            *np.indices(data.shape)) - data)
    res = optimize.leastsq(errorfunction, init, factor=factor,
                           full_output=True, maxfev=100000,
                           ftol=ftol, gtol=gtol, xtol=xtol)

    return res
