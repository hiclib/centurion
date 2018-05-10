# -*- coding: utf-8 -*-
from __future__ import print_function
import itertools
import numpy as np
from scipy import ndimage
from .optimization import fit_gaussian
from .externals import iced


def _detect_peaks(counts, lengths, sigma):
    """
    Detect peaks on sub_counts
    """
    peaks = []
    scores = []
    begin2, end2 = 0, 0
    for l2 in lengths:
        end2 += l2
        subcounts = ndimage.gaussian_filter(
            counts[:, begin2:end2], sigma)
        counts[:, begin2:end2] = subcounts
        begin2 = end2

    sum_trans = counts.sum(axis=1)
    der = sum_trans[1:] - sum_trans[:-1]

    if len(der) == 0:
        return peaks

    if der[0] < 0:
        peaks.append(0.5)
        scores.append(sum_trans[0])

    for i in range(1, len(der)):
        if np.sign(der[i]) < 0 and (np.sign(der[i - 1]) > 0) and \
           sum_trans[i] > np.median(sum_trans):
            peaks.append(i + 0.5)
            scores.append(- max(sum_trans[i - 1], sum_trans[i]) / 2)

    # Sort by score
    if len(scores) > 1:
        scores = np.array(scores)
        indx = scores.argsort()
        peaks = [peaks[i] for i in indx]
    return peaks


def find_centromeres_candidates(counts, lengths, n_candidate=3,
                                max_sigma=None, min_sigma=1, verbose=0,
                                return_sigma=False):
    """
    Find approximate centromeres positions by detecting peaks in the marginal
    of trans contact counts.

    Parameters
    ----------
    counts : ndarray n x n
        Contact count matrix.

    lengths : ndarray L
        Length of each chromosome.

    n_candidates : integer, optional, default: 3
        the min number of candidates.

    max_sigma : float, optional, default: None
        Maximum sigma to use. If None, is set to the chromosome length.

    min_sigma : float, optional, default : 1
        Minimum sigma to use.

    verbose : boolean, optional, default: False
        Verbose level

    return_sigma : boolean, optional, default: False
        If set to true, returns the L sigma in addition of the peaks

    Returns
    -------
    centromeres : a list of L list containing candidate centromeres for each
                  chromosome
    """
    if verbose:
        print("Searching for centromeres candidates")
    # FIXME if there is more than 3 peaks ?

    begin, end = 0, 0
    centromeres_call = []
    sigmas = []
    lencum = np.concatenate([[0], lengths.cumsum().astype(int)])

    counts = counts.copy()
    counts[:, lencum[:-1]] = 0
    counts[lencum[:-1]] = 0
    counts[:, lencum[1:] - 1] = 0
    counts[lencum[1:] - 1] = 0

    for length in lengths.astype(int):
        end += length
        sigma = 1
        max_sigma = length
        gf_counts = counts[begin:end]
        peaks = _detect_peaks(gf_counts.copy(), lengths, sigma)
        while (len(peaks) > n_candidate) and (sigma < max_sigma):
            sigma += 1
            gf_counts = counts[begin:end].copy()
            peaks = _detect_peaks(gf_counts.copy(), lengths, sigma)

        while len(peaks) < n_candidate and sigma >= min_sigma:
            sigma -= 0.5
            gf_counts = counts[begin:end].copy()
            peaks = _detect_peaks(gf_counts.copy(), lengths, sigma)

        # If there is no peaks, set the peak to the middle of the
        # centromere
        if len(peaks) == 0:
            peaks = [length / 2 + begin]
        else:
            peaks = [peak + begin for peak in peaks]
        centromeres_call.append(peaks)
        sigmas.append(sigma)
        begin = end
    if return_sigma:
        return centromeres_call, sigmas
    return centromeres_call


def filter_centromeres_candidates(counts, lengths, candidates, copy=True,
                                  sigma=4, verbose=0):
    """
    Filter centromeres candidates using a set of heuristics.

    Parameters
    ----------
    counts : ndarray n x n
        Contact count matrix.

    lengths : ndarray L
        Length of each chromosome.

    candidates : list of L liste
        list of candidates

    copy : boolean, optional, default: True
        whether to copy the contact counts matrix

    sigma : float, optional, default: 4
        initialization of the sigma parameters of the gaussians

    verbose : boolean, optional, default: False
        Verbose level

    Returns
    -------
    candidates : a list of L list containing the reduced candidate centromeres
                 for each chromosome.

    """
    if verbose:
        print("Filtering centromeres candidates")

    if copy:
        counts = counts.copy()
    baseline_candidates = [c[0] for c in candidates]
    mask = iced.utils.get_intra_mask(lengths)
    counts[mask] = 0
    parameters = np.concatenate([baseline_candidates,
                                 [counts.max() - np.median(counts),
                                  np.median(counts),
                                  sigma]])

    baseline_results, _, infodict, _, _ = fit_gaussian(
        counts, parameters, lengths, counts=True,
        factor=10, xtol=0.01)
    baseline = ((counts.flatten() - infodict["fvec"]) ** 2).sum()

    all_candidates = [[i, g] for i, cent in enumerate(candidates)
                      for j, g in enumerate(cent) if j != 0]

    kept_candidates = [[cent] for cent in baseline_results[:len(lengths)]]

    for num, candidate in all_candidates:
        parameters = baseline_results.copy()
        parameters[num] = candidate
        results, _, infodict, _, _ = fit_gaussian(
            counts, parameters, lengths, counts=True,
            factor=10, xtol=0.5)

        obj_value = ((counts.flatten() - infodict["fvec"]) ** 2).sum()
        if obj_value <= baseline:
            kept_candidates[num].append(candidate)
    return kept_candidates


def refine_centromeres(counts, lengths, candidate, sigma=4, verbose=0):
    """
    A single run of the optimization, assuming all preparation is complete
    """
    parameters = np.concatenate(
        [candidate,
         [counts.max() - np.median(counts),
          np.median(counts),
          sigma]])
    results, cov_x, infodict, mesg, suc = fit_gaussian(
        counts, parameters, lengths, counts=True,
        factor=10)
    fval = ((counts.flatten() - infodict["fvec"]) ** 2).sum()
    best_results = results[:len(lengths)] + 0.5
    return fval, best_results


def optimize_centromeres(counts, lengths, candidates, sigma=4, verbose=0,
                         njobs=1,
                         copy=True):
    """
    Perform the optimization steps of finding centromeres given a certain
    number of candidates.

    Parameters
    ----------
    counts : ndarray n x n
        Contact count matrix.

    lengths : ndarray L
        Length of each chromosome.

    candidates : list of L liste
        list of candidates

    copy : boolean, optional, default: True
        whether to copy the contact counts matrix

    sigma : float, optional, default: 4
        initialization of the sigma parameters of the gaussians


    """
    if verbose:
        print("Refining centromeres calls.")

    if copy:
        counts = counts.copy()

    mask = iced.utils.get_intra_mask(lengths)
    counts[mask] = 0

    all_candidates = itertools.product(
        *[candidate for candidate in candidates])

    if verbose:
        print("%d candidates" % np.prod([len(c) for c in candidates]))
        print

    fval_min = None
    best_results = None
    for i, c in enumerate(all_candidates):
        if verbose:
            print("Computing %d / %d candidates" % (
                i + 1, np.prod([len(j) for j in candidates])))
        fval, results = refine_centromeres(
            counts, lengths, c,
            sigma=sigma, verbose=verbose)
        if fval < fval_min or fval_min is None:
            fval_min = fval
            best_results = results

    return best_results
