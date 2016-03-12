#!/usr/bin/env python
from __future__ import print_function, division

import os
import os.path
import warnings

import yaml

import numpy as np
import matplotlib.pyplot as plt

from astropy.utils.compat import argparse
import astropy.table
import astropy.constants
import astropy.units as u

import desibest.utility


CLIGHT_KM_S = astropy.constants.c.to(u.km / u.s).value


def summarize(args=None):
    """Generate summary plots specified by a YAML config file.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action = 'store_true',
        help='provide verbose output on progress')
    parser.add_argument('config', type=str, metavar='YAML',
        help='Name of YAML configuration file to use')
    args = parser.parse_args(args)

    with open(args.config) as f:
        config = yaml.safe_load(f)

    title = config['title']
    if args.verbose:
        print('Generating summary plots for "{0}"'.format(title))
    nrows = config['rows']
    ncols = config['cols']
    base = config['base']
    find_file = lambda name: os.path.join(base, name)

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
        truth = astropy.table.Table.read(truth_name, hdu='_TRUTH')
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
                    zbest_dict[fitter_name][class_name], hdu='ZBEST')
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

            # Plot the performance of each fitter.
            for i in range(nrows * ncols):
                row = i // ncols
                col = i % ncols
                axis = axes[row][col]

                if i < len(fitter_names):
                    fitter_name = fitter_names[i]
                    ok = ok_dict[fitter_name]
                    dv = dv_dict[fitter_name]
                    bad = dv[ok] > catastrophic
                    if dv is not None:
                        lhs, rhs = desibest.utility.plot_slices(
                            x=x, y=dv, ok=ok, bad=bad, x_lo=node['min'],
                            x_hi=node['max'], num_slices=node['n'],
                            y_cut=max_dv, axis=axis)
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
                        # Why is -2 required here when -1 works below???
                        plt.setp([rhs.get_yticklabels()[-2:]], visible=False)
                    rhs.set_ylabel('zwarn, catastrophic fit fraction')

                if row < nrows - 1:
                    plt.setp([axis.get_xticklabels()], visible=False)
                else:
                    # Hide the last x-axis label except in the bottom right.
                    if col < ncols - 1:
                        plt.setp([axis.get_xticklabels()[-1:]], visible=False)
                    axis.set_xlabel('{0} {1}'.format(class_name, node['label']))

        figure.subplots_adjust(
            left=0.08, bottom=0.07, right=0.92, top=0.95, hspace=0., wspace=0.)
        plt.show()


if __name__ == '__main__':
    summarize()
