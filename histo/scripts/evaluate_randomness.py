#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from plot.plotting_utilities import *
import os
import sys
import pandas as pd
import argparse

plt.rcParams['font.family'] = 'DejaVu Sans Mono'
plt.rcParams['text.usetex'] = True


this_directory = os.path.dirname(os.path.realpath(__file__)) + "/"
this_filename = sys.argv[0].split('/')[-1]


res_dir = this_directory + '/../results/raw/randomness_test/'
images_dir = this_directory + '/../results/plots/randomness_test/'

if not os.path.exists(images_dir):
    os.makedirs(images_dir)

if not os.path.exists(res_dir):
    os.makedirs(res_dir)

parser = argparse.ArgumentParser(description='Plot the scalability of multiple implementations.',
                                 usage='')

parser.add_argument('-s', '--show', action='store_true',
                    help='Show the plots.')

parser.add_argument('-b', '--bins', type=int, default=1000,
                    help='Number of histogram bins.')
args = parser.parse_args()


plotconf = {
    'data_files': {
        'legacy std': [res_dir + 'n_p100M/synch_rad1_numbers.txt'],
        'legacy boost': [res_dir + 'n_p100M/synch_rad2_numbers.txt'],
        'std par': [
            res_dir + 'n_p100M/synch_rad4_0_numbers.txt',
            res_dir + 'n_p100M/synch_rad4_1_numbers.txt',
            res_dir + 'n_p100M/synch_rad4_2_numbers.txt',
            res_dir + 'n_p100M/synch_rad4_3_numbers.txt', ],
        'boost par': [
            res_dir + 'n_p100M/synch_rad6_0_numbers.txt',
            res_dir + 'n_p100M/synch_rad6_1_numbers.txt',
            res_dir + 'n_p100M/synch_rad6_2_numbers.txt',
            res_dir + 'n_p100M/synch_rad6_3_numbers.txt', ]
    },
    'hatches': ['', '', 'xx'],
    'markers': ['x', 'o', '^'],
    'xlabel': 'Bins',
    'ylabel': 'Density',
    # 'figsize': [5, 3],
    # 'xlim': [-6, 6],
    # 'ylim': [-0.01, 0.42],
    # 'outfiles': [images_dir+'randomness_test_100M_1.png',
    #              images_dir+'randomness_test_100M_1.pdf'],

    'figsize': [2.5, 2.5],
    'xlim': [-1, 1],
    'ylim': [0.3, 0.41],
    'outfiles': [images_dir+'randomness_test_100M_2.png',
                 images_dir+'randomness_test_100M_2.pdf'],

    'linestyles': [(0, (1, 3)), (1, (1, 3)), (2, (1, 3)), (3, (1, 3))],
    # 'colors': ['xkcd:red', 'xkcd:green', 'xkcd:blue'],
    # 'labels': ['legacy std', 'legacy boost',
    #            'std, thr0', 'std, thr1', 'std, thr2', 'std, thr3',
    #            'boost, thr0', 'boost, thr1', 'boost, thr2', 'boost, thr3',
    #            ],
    'title': {
        # 's': '{}'.format(case.upper()),
        'fontsize': 10,
        # 'y': .85,
        # 'x': 0.55,
        # 'fontweight': 'bold',
    },
    'annotate': {
        'fontsize': 9,
        'textcoords': 'data',
        'va': 'bottom',
        'ha': 'center'
    },
    'ticks': {'fontsize': 10},
    'fontsize': 10,
    'legend': {
        'loc': 'upper left', 'ncol': 1, 'handlelength': 2, 'fancybox': False,
        'framealpha': .9, 'fontsize': 10, 'labelspacing': 0.5, 'borderpad': 0.5,
        'handletextpad': 0.5, 'borderaxespad': 0.2, 'columnspacing': 0.8,
        # 'bbox_to_anchor': (0., 0.85)
    },
    'subplots_adjust': {
        'wspace': 0.0, 'hspace': 0.1, 'top': 0.93
    },
    'tick_params': {
        'pad': 1, 'top': 0, 'bottom': 1, 'left': 1,
        'direction': 'out', 'length': 3, 'width': 1,
    },
}


def import_results(files):
    d = np.array([])
    for file in files:
        temp = pd.read_csv(file, delimiter='\n', dtype=float)
        d = np.concatenate((d, temp.to_numpy().reshape((-1,))))
    return d


if __name__ == '__main__':

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=plotconf['figsize'])
    i = 0
    for label, data_files in plotconf['data_files'].items():
        x = import_results(data_files)
        hist, edges = np.histogram(x, bins=args.bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        plt.plot(centers, hist, label=label, lw=2, ls=plotconf['linestyles'][i])
        i+= 1
    plt.xlim(plotconf['xlim'])
    plt.ylim(plotconf['ylim'])
    plt.xticks(**plotconf['ticks'])
    plt.yticks(**plotconf['ticks'])
    ax.tick_params(**plotconf['tick_params'])
    if plotconf['xlim'][0] != -1:
        ax.legend(**plotconf['legend'])
        plt.ylabel(plotconf['ylabel'], fontsize=plotconf['fontsize'])
        plt.xlabel(plotconf['xlabel'], fontsize=plotconf['fontsize'])

    plt.tight_layout()
    plt.subplots_adjust(**plotconf['subplots_adjust'])
    for file in plotconf['outfiles']:
        save_and_crop(fig, file, dpi=600, bbox_inches='tight')
    if args.show:
        plt.show()
    plt.close()
