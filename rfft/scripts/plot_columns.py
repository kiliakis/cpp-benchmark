#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np

res_dir = '/afs/cern.ch/work/k/kiliakis/results/csv/'
images_dir = '/afs/cern.ch/work/k/kiliakis/results/images/'

csv_files = [res_dir + 'cpp/fft/v1/fft-v1-fft.csv',
             res_dir + 'cpp/fft/v1/fft-v1-ifft.csv',
             res_dir + 'cpp/fft/v1/fft-v1-rfft.csv',
             res_dir + 'cpp/fft/v1/fft-v1-irfft.csv']

x_names = ['n_points'] * 4
y_names = ['Throughput'] * 4
err_name = 'Stdev'
line_names = ['fft', 'ifft', 'rfft', 'irfft']
image_name = images_dir+'fft-v1.png'
x_label = 'Points'
y_label = 'Throughput [MB/s]'
title = 'FFTs Throughput v1'

x_lims = []
y_lims = [0, 1400]

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
    d = np.loadtxt(file, dtype='str')
    return d.tolist()


def plot(x, y, label, xerr=None, yerr=None):

    plt.grid(True, which='major', alpha=1)
    plt.grid(True, which='minor', alpha=1)
    plt.minorticks_on()
    plt.errorbar(x, y, xerr=xerr, yerr=yerr, marker='o', linewidth='2', label=label)


if __name__ == '__main__':

    plt.figure(figsize=(10, 6))
    if(x_lims):
        plt.xlim(x_lims)
    if(y_lims):
        plt.ylim(y_lims)
    #ax = plt.gca()
    #ax.set_xlim(x_lims)
    #ax.set_ylim(y_lims)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tick_params(labelright=True)
    plt.title(title)
    plt.xscale('log')

    for file, x_name, y_name, line_name in zip(csv_files, x_names, y_names, line_names):
        l = import_results(file)
        header = l[0]
        array = np.array(l[1:])
        c = header.index(x_name)
        x = array[:, c].astype(float)
        c = header.index(y_name)
        y = array[:, c].astype(float)

        if err_name:
            c = header.index(err_name)
            yerr = array[:, c].astype(float)
            plot(x, y, line_name, yerr=yerr)
        else:
            plot(x, y, line_name)

        if annotate_min_flag:
            annotate_min(x, y)
        if annotate_max_flag:
            annotate_max(x, y)

    plt.legend(loc='best', fancybox=True, framealpha=0.5)
    plt.savefig(image_name, bbox_inches='tight')
    plt.tight_layout()
    plt.close()
