#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os

this_directory = os.path.dirname(os.path.realpath(__file__)) + "/"

res_dir = this_directory + '../results/csvfiles/'
images_dir = this_directory + '../results/plots/'

csv_file = res_dir + 'v2.csv'

image_name = images_dir+'v2.pdf'
x_label = 'Number of Points [1e6]'
y_label = 'Throughput [MPoints/s]'
title = 'Drift benchmark'
dict_keys = ['alpha', 'n_threads']
x_key = 'n_points'
y_key = 'throughput(mp/sec)'

x_lims = []
y_lims = [0, 1600]

annotate_min_flag = False
annotate_max_flag = False


def annotate(ax, A, B):
    for xy in zip(A, B):
        ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data')


def annotate_min(A, B):
    y = min(B)
    i = B.index(y)
    plt.subplot().annotate('%.2e' % float(y), xy=(A[i], y), size='large')


def annotate_max(A, B):
    y = max(B)
    i = B.index(y)
    plt.subplot().annotate('%.2e' % float(y), xy=(A[i], y), size='large')


def import_results(file):
    d = np.genfromtxt(file, dtype='str', delimiter='\t')
    return d.tolist()


def plot(x, y, label, xerr=None, yerr=None):

    plt.grid(True, which='major', alpha=1)
    plt.grid(True, which='minor', alpha=1)
    plt.minorticks_on()
    plt.errorbar(x, y, xerr=xerr, yerr=yerr,
                 marker='o', linewidth='2', label=label)


if __name__ == '__main__':
    data = import_results(csv_file)
    header, data = list(data[0]), data[1:]
    dic = {}
    for r in data:
        if r[header.index(dict_keys[0])] not in dic:
            dic[r[header.index(dict_keys[0])]] = {}
        if r[header.index(dict_keys[1])] not in dic[r[header.index(dict_keys[0])]]:
            dic[r[header.index(dict_keys[0])]][r[header.index(dict_keys[1])]] = {'x': [], 'y': []}
        dic[r[header.index(dict_keys[0])]][r[header.index(dict_keys[1])]]['x'].append(int(r[header.index(x_key)]))
        dic[r[header.index(dict_keys[0])]][r[header.index(dict_keys[1])]]['y'].append(float(r[header.index(y_key)]))

    plt.figure(figsize=(10, 6))

    # if(x_lims):
    #     plt.xlim(x_lims)
    # if(y_lims):
    #     plt.ylim(y_lims)
    #ax = plt.gca()
    # ax.set_xlim(x_lims)
    # ax.set_ylim(y_lims)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    # plt.tick_params(labelright=True)
    plt.title(title)
    plt.xscale('log')
    plt.grid(axis='both')
    config = {
        '0': {'marker': 'o', 'ls': '-'},
        '1': {'marker': 'o', 'ls': '-'},
        '2': {'marker': 'o', 'ls': '-'},
        'legacy': {'marker': 'o', 'ls': '-'},
        'new': {'marker': 's', 'ls': '--'}
    }

    for bench, v1 in dic.items():
        colors = ['blue', 'orange', 'green', 'red']
        for thr, xy in v1.items():
            if thr in ['4', '56']:
                continue
            x = np.array(xy['x'])//1000000
            y = np.array(xy['y'])
            idx = np.argsort(x)
            x, y = x[idx], y[idx]
            plt.plot(x, y, lw=2, color=colors.pop(),
                     marker=config[bench]['marker'],
                     ls=config[bench]['ls'],
                     label='{}-{}thr'.format(bench, thr))
    plt.xticks(x, x)

    # for file, x_name, y_name, line_name in zip(csv_files, x_names, y_names, line_names):
    #     l = import_results(file)
    #     header = l[0]
    #     array = np.array(l[1:])
    #     c = header.index(x_name)
    #     x = array[:, c].astype(float)
    #     c = header.index(y_name)
    #     y = array[:, c].astype(float)

    #     if err_name:
    #         c = header.index(err_name)
    #         yerr = array[:, c].astype(float)
    #         plot(x, y, line_name, yerr=yerr)
    #     else:
    #         plot(x, y, line_name)

    #     if annotate_min_flag:
    #         annotate_min(x, y)
    #     if annotate_max_flag:
    #         annotate_max(x, y)

    plt.tight_layout()
    plt.legend(loc='upper right', fancybox=True, framealpha=0.5,
               ncol=2)
    plt.savefig(image_name, bbox_inches='tight', dpi=600)
    plt.show()
    plt.close()
