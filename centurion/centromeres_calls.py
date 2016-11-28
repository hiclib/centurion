import numpy as np
from . import prelocalization_
from . import utils
from .externals import iced


def centromeres_calls(counts, lengths, resolution=40000, init=None,
                      n_candidate=3,
                      sigma=None, verbose=0, filter_candidates=True,
                      normalize=True):
    """
    Calls centromeres

    Parameters
    """
    if sigma is None:
        sigma = 80000 / resolution

    coef = 40000 / resolution

    if init is None:
        counts_40kb, lengths_40kb = utils.downsample_resolution(
            counts, lengths, coefficient=coef)

        if normalize:
            counts_40kb = iced.normalization.ICE_normalization(counts_40kb)

        centromeres_calls_40kb = centromeres_calls_(
            counts_40kb, lengths_40kb, sigma=sigma, verbose=verbose,
            filter_candidates=filter_candidates)

        candidates = []
        for i, candidate in enumerate(centromeres_calls_40kb):
            candidates.append(
                [(candidate - lengths_40kb[:i].sum()) * coef +
                 lengths[:i].sum()])

    else:
        init /= resolution
        if len(init) > 1:
            init[1:] += lengths[:-1].cumsum()
        candidates = [[i] for i in init]

    if normalize:
        counts = iced.normalization.ICE_normalization(counts)

    mask = iced.utils.get_intra_mask(lengths)
    counts[mask] = 0
    if coef != 1:
        centromeres_calls = prelocalization_.optimize_centromeres(
            counts, lengths, candidates,
            sigma=1, verbose=verbose)
    else:
        centromeres_calls = centromeres_calls_40kb
    centromeres_calls[1:] -= lengths[:-1].cumsum()
    centromeres_calls *= resolution
    return np.round(centromeres_calls).astype(int)


def centromeres_calls_(counts, lengths, sigma=4, init=None,
                       n_candidate=2, verbose=2, copy=True,
                       filter_candidates=False, max_trials=30):
    """
    """
    if copy:
        counts = counts.copy()
    mask = iced.utils.get_intra_mask(lengths)
    counts[mask] = 0

    candidate_centromeres = prelocalization_.find_centromeres_candidates(
        counts, lengths, n_candidate=n_candidate, verbose=verbose)
    potential_candidates = np.prod([len(i) for i in candidate_centromeres])
    if filter_candidates and potential_candidates > max_trials:
        refined_candidates = prelocalization_.filter_centromeres_candidates(
            counts, lengths, candidate_centromeres, sigma=sigma,
            verbose=verbose)
    else:
        refined_candidates = candidate_centromeres
    centromeres_calls = prelocalization_.optimize_centromeres(
        counts, lengths, refined_candidates, sigma=sigma, verbose=verbose)
    return centromeres_calls
