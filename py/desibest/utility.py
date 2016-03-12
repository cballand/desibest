"""Utility functions.
"""
from __future__ import print_function, division

import numpy as np


def plot_slices(x, y, x_lo, x_hi, y_cut, num_slices=5, min_count=100,
                axis=None):
    """Scatter plot with 68, 95 percentiles superimposed in slices.

    Requires that the matplotlib package is installed.

    Parameters
    ----------
    x : array
        X-coordinates to scatter plot.  Points outside [x_lo, x_hi] are
        not displayed.
    y : array
        Y-coordinates to scatter plot.  Y values are assumed to be roughly
        symmetric about zero.
    x_lo : float
        Minimum value of x to plot.
    x_hi : float
        Maximum value of x to plot.
    y_cut : float
        The target maximum value of |y|.  A dashed line at this value is
        added to the plot, and the vertical axis is clipped at
        |y| = 1.25 * y_cut (but values outside this range are included in
        the percentile statistics).
    num_slices : int
        Number of equally spaced slices to divide the interval [x_lo, x_hi]
        into.
    min_count : int
        Do not use slices with fewer points for superimposed percentile
        statistics.
    axis : matplotlib axis object or None
        Uses the current axis if this is None.
    """
    import matplotlib.pyplot as plt

    if axis is None:
        axis = plt.gca()

    # Assume that y is symmetric about zero.
    axis.set_xlim(x_lo, x_hi)
    axis.set_ylim(-1.25 * y_cut, +1.25 * y_cut)
    x_bins = np.linspace(x_lo, x_hi, num_slices + 1)
    x_i = np.digitize(x, x_bins) - 1
    limits = []
    counts = []
    for s in xrange(num_slices):
        y_slice = y[x_i == s]
        counts.append(len(y_slice))
        if counts[-1] > 0:
            limits.append(np.percentile(y_slice, (2.5, 16, 50, 84, 97.5)))
        else:
            limits.append((0., 0., 0., 0., 0.))
    limits = np.array(limits)
    counts = np.array(counts)

    # Plot scatter of all fits.
    axis.scatter(x, y, s=15, marker='.', lw=0, color='blue', alpha=0.5)

    # Plot quantiles in slices with enough fits.
    stepify = lambda y: np.vstack([y, y]).transpose().flatten()
    y_m2 = stepify(limits[:, 0])
    y_m1 = stepify(limits[:, 1])
    y_med = stepify(limits[:, 2])
    y_p1 = stepify(limits[:, 3])
    y_p2 = stepify(limits[:, 4])
    xstack = stepify(x_bins)[1:-1]
    for i in xrange(num_slices):
        s = slice(2 * i, 2 * i + 2)
        if counts[i] >= min_count:
            axis.fill_between(
                xstack[s], y_m2[s], y_p2[s], alpha=0.15, color='red')
            axis.fill_between(
                xstack[s], y_m1[s], y_p1[s], alpha=0.25, color='red')
            axis.plot(xstack[s], y_med[s], 'r-', lw=2.)

    # Plot cut lines.
    axis.axhline(+y_cut, ls='-', color='k')
    axis.axhline(-y_cut, ls='-', color='k')
    axis.grid()
