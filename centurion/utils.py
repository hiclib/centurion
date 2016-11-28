import numpy as np


def downsample_resolution(counts, lengths, coefficient=2):
    """
    Downsamples the resolution of a matrix

    Parameters
    ----------
    counts : ndarray (N, N)
        contact counts matrix to downsample

    lengths : ndarray (L, )
        chromosomes lengths

    coef : int, optionnal, default: 2
        downsample resolution of the counts matrix by `coef`

    Returns
    -------
    target_counts, target_lengths : ndarray
    """
    if coefficient == 1:
        return counts, lengths
    # FIXME there is probably a better way to do this
    target_lengths = np.ceil(lengths.astype(float) / coefficient).astype(int)
    target_counts = np.zeros((target_lengths.sum(),
                              target_lengths.sum()))
    begin_i, end_i = 0, 0
    target_begin_i, target_end_i = 0, 0
    for i, length_i in enumerate(lengths):
        end_i += length_i
        target_end_i += target_lengths[i]
        begin_j, end_j = 0, 0
        target_begin_j, target_end_j = 0, 0
        for j, length_j in enumerate(lengths):
            end_j += length_j
            target_end_j += target_lengths[j]

            sub_counts = counts[begin_i:end_i, begin_j:end_j]
            sub_target_counts = target_counts[target_begin_i:target_end_i,
                                              target_begin_j:target_end_j]
            for start in range(coefficient):
                s = sub_counts[start::coefficient, start::coefficient]
                sub_target_counts[:s.shape[0], :s.shape[1]] += s

            begin_j = end_j
            target_begin_j = target_end_j
        begin_i = end_i
        target_begin_i = target_end_i
    return target_counts, target_lengths
