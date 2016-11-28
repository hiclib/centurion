import numpy as np
from centurion.externals import iced
from centurion import centromeres_calls
import matplotlib.pyplot as plt
from matplotlib import colors

counts, lengths = iced.datasets.load_sample_yeast()
centromeres = centromeres_calls.centromeres_calls(
    counts, lengths,
    resolution=10000)

# Normalize the data for the sake of visualization
counts = iced.filter.filter_low_counts(counts, percentage=0.04)
counts = iced.normalization.ICE_normalization(counts)

# And remove the intra chromosomal for the sake of visualization
mask = iced.utils.get_intra_mask(lengths)
counts[mask] = np.nan

fig, ax = plt.subplots()
ax.matshow(counts, cmap="RdBu", norm=colors.LogNorm())
centro = centromeres / 10000
centro[1:] += lengths.cumsum()[:-1]
[ax.axhline(i, color="#000000", linestyle="--") for i in centro]
[ax.axvline(i, color="#000000", linestyle="--") for i in centro]
