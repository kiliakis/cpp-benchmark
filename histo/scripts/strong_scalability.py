import matplotlib.pyplot as plt
import numpy as np

import os
from matplotlib import cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from itertools import cycle
import matplotlib.ticker
import sys
from plot.plotting_utilities import *
import argparse

this_directory = os.path.dirname(os.path.realpath(__file__)) + "/"
this_filename = sys.argv[0].split('/')[-1]

parser = argparse.ArgumentParser(description='Plot the scalability of multiple implementations.',
                                 usage='')

# parser.add_argument('-c', '--cases', type=str, nargs='+',
#                     choices=['lhc', 'sps', 'ps', 'ex01'],
#                     help='The test-case to plot.')

# parser.add_argument('-k', '--keysuffix', type=str, default='strong',
#                     help='A key suffix to use.')

# parser.add_argument('-e', '--errorbars', action='store_true',
#                     help='Add errorbars.')

parser.add_argument('-s', '--show', action='store_true',
                    help='Show the plots.')

args = parser.parse_args()

project_dir = this_directory + '../'
res_dir = project_dir + 'results/'
images_dir = res_dir + 'plots/'

if not os.path.exists(images_dir):
    os.makedirs(images_dir)

gconfig = {
    'x_name': 'n_threads',
    'x_to_keep': [1, 2, 8, 14, 28],
    'y_name': 'time',
    'y_err_name': 'stdev',
    'fontname': 'DejaVu Sans Mono',

    # 'ylim': [0, 28],
}

label_d = {
    'smooth_histo1': 'smooth-orig',
    'smooth_histo4': 'smooth-opt',
    # 'histo': '',
    # 'synch_rad4': 'std, par, single',
    # 'synch_rad5': 'boost, par, sep',
    # 'synch_rad6': 'boost, par, single',
    # 'synch_rad7': 'std, par, single, vec',
    # 'synch_rad8': 'boost, par, single, vec',
    # 'synch_rad9': 'boost, par, single, vec, opt',
}

lconfig = {
    'figures': {
        'strong': {
            'files': [
                '{}/csv/perftest.csv',
            ],
            'lines': {
                'exe': ['smooth_histo1', 'smooth_histo4'],
                'n_points': ['1000000', '4000000', '16000000'],
            }
        },
    },
}

tabconf = {
    'exe': ['smooth_histo1', 'smooth_histo4'],
    'n_points': '16000000',
    'outfile': '{}/table2.csv'
}

plotconf = {
    'annotate_exe': ['smooth_histo1', 'smooth_histo4'],
    'exe': ['smooth_histo1', 'smooth_histo4'],
    # 'n_points': ['1000000', '2000000', '4000000', '8000000', '16000000', '32000000'],
    'n_points': ['16000000'],
    'outfiles': ['{}/n_p{}scalability2.png',
                 '{}/n_p{}scalability2.pdf'],
    'ylim': [0, 36],
    # 'xlim': [, 36],
    'yticks': [4, 8, 12, 16, 20, 24, 28, 32],
    'figsize': [5, 3],
    'hatches': ['', '', 'xx'],
    'markers': ['x', 'o', '^'],
    'colors': ['xkcd:red', 'xkcd:green', 'xkcd:blue'],
    'xlabel': 'Cores',
    'ylabel': 'Throughput',
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
        'loc': 'upper left', 'ncol': 1, 'handlelength': 1.2, 'fancybox': False,
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

plt.rcParams['font.family'] = gconfig['fontname']
plt.rcParams['text.usetex'] = True

if __name__ == '__main__':
    plots_dir = {}
    for title, figconf in lconfig['figures'].items():
        for file in figconf['files']:
            # file = file.format(res_dir, case.upper())
            # print(file)
            data = np.genfromtxt(file.format(res_dir),
                                 delimiter='\t', dtype=str)
            header, data = list(data[0]), data[1:]
            temp = get_plots(header, data, figconf['lines'],
                             exclude=figconf.get('exclude', []),
                             prefix=True)
            for key in temp.keys():
                plots_dir['_{}_'.format(key)] = temp[key].copy()

    exe_names = []
    err_l = []
    times_l = []
    for key, vals in plots_dir.items():
        exe = key.split('exe')[1].split('_n_points')[0]
        n_points = key.split('n_points')[1].split('_')[0]
        if exe in tabconf['exe'] and n_points == tabconf['n_points']:
            # Here I want to keep only with one thread
            threads = get_values(vals, header, gconfig['x_name'])
            idx = threads.tolist().index(1)
            err = get_values(vals, header, gconfig['y_err_name'])[idx]
            time = get_values(vals, header, gconfig['y_name'])[idx]
            # Then I want to print a nice table with exe name, points, time, stdev, speedup
            exe_names.append(exe)
            err_l.append(err)
            times_l.append(time)

    err_l = np.array(err_l, dtype=float)
    times_l = np.array(times_l, dtype=float)
    err_l = 100 * err_l / times_l
    sortidx = np.argsort(exe_names)
    exe_names = np.array(exe_names)[sortidx]
    err_l = err_l[sortidx]
    times_l = times_l[sortidx]
    speedups = times_l[0] / times_l
    str = '{:<12}\t{:<6}\t{:<4}\t{:<4}\n'.format(
        'exe', 'time', 'stdev', 'speedup')
    for exe, t, e, s in zip(exe_names, times_l, err_l, speedups):
        str += '{:<12}\t{:<6.5}\t{:<4.2}\t{:<4.2}\n'.format(exe, t, e, s)
    print(str)
    print(str, file=open(tabconf['outfile'].format(images_dir), 'w'))
    # sys.exit()

    fig, ax_arr = plt.subplots(ncols=1, nrows=1,
                               sharex=True, sharey=True,
                               figsize=plotconf['figsize'])
    for key, values in plots_dir.items():
        exe = key.split('exe')[1].split('_n_points')[0]
        n_points = key.split('n_points')[1].split('_')[0]
       

        label = label_d[exe]
        # label = exe

        x = get_values(values, header, gconfig['x_name'])
        y = get_values(values, header, gconfig['y_name'])
        yerr = get_values(values, header, gconfig['y_err_name'])
        if len(x) == 1:
            x = np.array([x[0]] * len(gconfig['x_to_keep']))
            y = np.array([y[0]] * len(gconfig['x_to_keep']))
            yerr = np.array([yerr[0]] * len(gconfig['x_to_keep']))
        yerr = yerr / y
        y = int(n_points) * 100 / y
        yerr = y * yerr

        # speedup = y / yref
        # yerr = yerr / yref

        plt.errorbar(np.arange(len(x)), y,
                     label=label,
                     # marker=plotconf['markers'][idx],
                     # color=plotconf['colors'][idx],
                     yerr=yerr,
                     linewidth=1,
                     # linestyle='--',
                     capsize=2,
                     zorder=1)
        if exe in plotconf['annotate_exe']:
            for xi, yi in zip(np.arange(len(x)), y):
                plt.gca().annotate(r'\bfseries {:.1f}x'.format(yi/y[0]), xy=(xi, yi+0.4),
                            **plotconf['annotate'],
                            zorder=2)
    
    plt.grid(axis='y', zorder=0, alpha=0.75)
    plt.gca().set_axisbelow(True)

    # plt.grid(True, which='both', axis='y', alpha=0.9)
    # plt.grid(False, which='major', axis='x')
    # plt.title('{} Particles'.format(
    #     human_format(int(n_points))), **plotconf['title'])

    plt.xlabel(plotconf['xlabel'], labelpad=3,
               fontsize=plotconf['fontsize'])
    plt.ylabel(plotconf['ylabel'], labelpad=3,
               fontsize=plotconf['fontsize'])
    plt.xticks(np.arange(len(x)), np.array(
        x, dtype=int), **plotconf['ticks'])

    plt.gca().tick_params(**plotconf['tick_params'])

    plt.gca().legend(**plotconf['legend'])

    plt.xticks(**plotconf['ticks'])
    plt.yticks(**plotconf['ticks'])
    # plt.yticks(plotconf['yticks'], plotconf['yticks'], **plotconf['ticks'])

    plt.tight_layout()
    plt.subplots_adjust(**plotconf['subplots_adjust'])
    for file in plotconf['outfiles']:
        file = file.format(images_dir, human_format(int(n_points)))
        save_and_crop(fig, file, dpi=300, bbox_inches='tight')
    if args.show:
        plt.show()
    plt.close()
