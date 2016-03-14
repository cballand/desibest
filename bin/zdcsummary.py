#!/usr/bin/env python
from __future__ import print_function, division

import os
import os.path
import warnings
import sys

import yaml

import numpy as np

import matplotlib as mpl

from astropy.utils.compat import argparse
import astropy.table
import astropy.constants
import astropy.units as u


CLIGHT_KM_S = astropy.constants.c.to(u.km / u.s).value


def plot_slices(x, y, ok, bad, x_lo, x_hi, y_cut, num_slices=5, min_count=100,
                axis=None):
    """Scatter plot with 68, 95 percentiles superimposed in slices.

    Requires that the matplotlib package is installed.

    Parameters
    ----------
    x : array of float
        X-coordinates to scatter plot.  Points outside [x_lo, x_hi] are
        not displayed.
    y : array of float
        Y-coordinates to scatter plot.  Y values are assumed to be roughly
        symmetric about zero.
    ok : array of bool
        Array of booleans that identify which fits are considered good.
    bad : array of bool
        Array of booleans that identify which fits have failed catastrophically.
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

    x_bins = np.linspace(x_lo, x_hi, num_slices + 1)
    x_i = np.digitize(x, x_bins) - 1
    limits = []
    counts = []
    for s in xrange(num_slices):
        # Calculate percentile statistics for ok fits.
        y_slice = y[ok & (x_i == s)]
        counts.append(len(y_slice))
        if counts[-1] > 0:
            limits.append(np.percentile(y_slice, (2.5, 16, 50, 84, 97.5)))
        else:
            limits.append((0., 0., 0., 0., 0.))
    limits = np.array(limits)
    counts = np.array(counts)

    # Plot scatter of all fits.
    axis.scatter(x[ok], y[ok], s=15, marker='.', lw=0, color='b', alpha=0.5)
    axis.scatter(x[~ok], y[~ok], s=15, marker='x', lw=0, color='k', alpha=0.5)

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
    axis.axhline(+y_cut, ls=':', color='k')
    axis.axhline(0., ls='-', color='k')
    axis.axhline(-y_cut, ls=':', color='k')

    # Plot histograms of of not ok and catastrophic fits.
    rhs = axis.twinx()

    weights = np.ones_like(x[bad]) / len(x[ok])
    if len(weights) > 0:
        rhs.hist(
            x[bad], range=(x_lo, x_hi), bins=num_slices, histtype='step',
            weights=weights, color='k')

    weights = np.ones_like(x[~ok]) / len(x)
    if len(weights) > 0:
        rhs.hist(
            x[~ok], range=(x_lo, x_hi), bins=num_slices, histtype='step',
            weights=weights, color='k', ls='dashed')

    axis.set_ylim(-1.25 * y_cut, +1.25 * y_cut)
    axis.set_xlim(x_lo, x_hi)

    return axis, rhs


def summarize(args=None):
    """Generate summary plots specified by a YAML config file.

    See the example YAML files for details on how to configure the summary
    plots.

    This is an entry point for a command-line program.  The program runs in
    batch mode, saving all plots to a single multi-page pdf file, when the
    --output option is specified.  Otherwise, plots are displayed interactively
    and the plot window must be closed after page is displayed in order to
    advance to the next page.

    Parameter
    ---------
    args : list or None
        Command-line arguments to use.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true',
        help='provide verbose output on progress')
    parser.add_argument('-o', '--output', type=str, default=None,
        help='save plots to the specified pdf file instead of displaying them')
    parser.add_argument('config', type=str, metavar='YAML',
        help='Name of YAML configuration file to use')
    args = parser.parse_args(args)

    with open(args.config) as f:
        config = yaml.safe_load(f)

    if args.output:
        _, ext = os.path.splitext(args.output)
        if ext not in ('', '.pdf'):
            print('Can only save plots to a .pdf file')
            return -1
        if ext == '':
            args.output += '.pdf'
        mpl.use('pdf')
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_output = PdfPages(args.output)
    else:
        pdf_output = None

    # Defer this import until after we have selected the backend to use.
    import matplotlib.pyplot as plt

    title = config['title']
    if args.verbose:
        print('Generating summary plots for "{0}"'.format(title))
    nrows = config['rows']
    ncols = config['cols']
    base = config['base']

    # Define a helper function to find files used below.
    def find_file(name):
        full_name = os.path.join(base, name)
        if not os.path.isfile(full_name):
            print('No such file: {0}'.format(full_name))
            sys.exit(-1)
        return full_name

    zbest_dict = config['zbest']
    fitter_names = zbest_dict.keys()
    assert len(fitter_names) <= nrows * ncols

    # Silence expected matplotlib warnings.
    warnings.simplefilter('ignore', category=FutureWarning)
    warnings.simplefilter('ignore', category=UserWarning)

    # Loop over target classes.
    class_names = config['classes'].keys()
    for class_name in class_names:

        node = config['classes'][class_name]
        max_dv = node['max_dv']
        catastrophic = node['catastrophic']
        max_frac = node['max_frac']

        # Load the truth information.
        truth_name = node['truth']
        truth = astropy.table.Table.read(find_file(truth_name), hdu='_TRUTH')
        z_true = truth['TRUEZ'].data
        if args.verbose:
            print(
                'Loaded {0} columns of {1} truth from {2}.'
                .format(len(truth), class_name, truth_name))

        # Look for results from each fitter for this class.
        subplot_index = 0
        dv_dict = {}
        ok_dict = {}
        for fitter_name in fitter_names:
            if class_name in zbest_dict[fitter_name]:
                zbest = astropy.table.Table.read(
                    find_file(zbest_dict[fitter_name][class_name]), hdu='ZBEST')
                # Calculate velocity residuals.
                z_est = zbest['Z'].data
                ok = (zbest['ZWARN'].data == 0)
                ok_dict[fitter_name] = ok
                dv = CLIGHT_KM_S * (z_est - z_true) / (1 + z_true)
                dv_dict[fitter_name] = dv
                if args.verbose:
                    dv_quant = np.percentile(dv[ok], (2.5, 16, 50, 84, 97.5))
                    dv_err68 = (dv_quant[3] - dv_quant[1])/2.
                    dv_err95 = (dv_quant[4] - dv_quant[0])/2.
                    ok_pct = 1e2 * np.count_nonzero(ok) / len(ok)
                    print(
                        '{0} {1}: med {2:.1f} 68% {3:.1f} 95% {4:.1f} km/s'
                        ', ok {5:.1f}%'
                        .format(class_name, fitter_name, dv_quant[2],
                                dv_err68, dv_err95, ok_pct))
            else:
                ok_dict[fitter_name] = None
                dv_dict[fitter_name] = None

        # Loop over plots to make for this target class.
        for plot_var_name, node in node['plots'].iteritems():

            # Initialize a new page of plots.
            figure, axes = plt.subplots(
                nrows, ncols, figsize=(11, 8.5), facecolor='white',
                sharex=True, sharey=True)
            figure.suptitle(title)

            # Plot the truth distribution for this variable.
            x = truth[plot_var_name].data
            if args.verbose:
                print('Plotting {0} {1}'.format(class_name, plot_var_name))
            n, x_min, x_max = node['n'], node['min'], node['max']
            overlap = node['overlap']

            # Plot the performance of each fitter.
            for i in range(nrows * ncols):
                row = i // ncols
                col = i % ncols
                axis = axes[row][col]

                if i < len(fitter_names):
                    fitter_name = fitter_names[i]
                    ok = ok_dict[fitter_name]
                    dv = dv_dict[fitter_name]
                    if dv is not None:
                        bad = dv[ok] > catastrophic
                        lhs, rhs = desibest.utility.plot_slices(
                            x=x, y=dv, ok=ok, bad=bad, x_lo=x_min, x_hi=x_max,
                            num_slices=n, y_cut=max_dv, axis=axis)
                    # Add a label even if the fitter has no results.
                    xy = (0.5, 1.0)
                    coords = 'axes fraction'
                    axis.annotate(
                        fitter_name, xy=xy, xytext=xy, xycoords=coords,
                        textcoords=coords, horizontalalignment='center',
                        verticalalignment='top', size='large', weight='bold')

                else:
                    rhs = axis.twinx()

                rhs.set_ylim(0., max_frac)
                if col > 0:
                    plt.setp([axis.get_yticklabels()], visible=False)
                else:
                    axis.set_ylabel('Redshift fit residual $\Delta v$ [km/s]')

                if col < ncols - 1:
                    plt.setp([rhs.get_yticklabels()], visible=False)
                else:
                    # Hide the last y-axis label except on the first row.
                    if row > 0:
                        # Why is -2 required here??
                        plt.setp([rhs.get_yticklabels()[-2:]], visible=False)
                    rhs.set_ylabel('zwarn, catastrophic fit fraction')

                if row < nrows - 1:
                    plt.setp([axis.get_xticklabels()], visible=False)
                else:
                    axis.set_xlabel('{0} {1}'.format(class_name, node['label']))
                    # Hide overlapping x-axis labels except in the bottom right.
                    if overlap and (col < ncols - 1):
                        plt.setp(
                            [axis.get_xticklabels()[-overlap:]], visible=False)

            figure.subplots_adjust(
                left=0.08, bottom=0.07, right=0.92, top=0.95,
                hspace=0., wspace=0.)

            if pdf_output:
                pdf_output.savefig()
                plt.close()
            else:
                print('Close the plot window to continue...')
                sys.stdout.flush()
                plt.show()

    if pdf_output:
        meta = pdf_output.infodict()
        meta['Title'] = title
        pdf_output.close()
        if args.verbose:
            print('Saved plots to {0}'.format(args.output))


if __name__ == '__main__':
    summarize()
