#!/usr/bin/python
#import os
#import csv
import numpy as np

from prettytable import PrettyTable

input_file = '/afs/cern.ch/work/k/kiliakis/results/csv/cpp/' +\
            'fft/v1/fft-v1.csv'

output_file_prefix = '/afs/cern.ch/work/k/kiliakis/results/csv/cpp/' +\
            'fft/v1/fft-v1'
column_names = ['fft_type']
postfixes = ['-']
#x_names = ['n_threads', 'n_threads']
#y_names = ['turn_mean', 'turn_mean']
#line_names = ['cpp 1M partciles', 'python 1M partciles']
#image_name = 'image.png'
#x_label = 'Threads'
#y_label = 'Average Turn Time (sec)'
#title = 'Python - C++ Run Time Comparison'

header = []

def group_by(list, name, prefix):
    global header
    header = list[0]
    #array = np.array(list[1:])
    j = header.index(name)
    d = {}
    for r in list[1:]:
        key = prefix + r[j]
        if key in d:
            d[key].append(r)
        else:
            d[key] = [r]
    return d
    #column = array[:, j]
    #unique_values = np.unique(column)


def dump_to_files(postfix, records):
    t = PrettyTable(header)
    t.align = 'l'
    t.border = False
    for r in records:
        t.add_row(r)
    with open(output_file_prefix+postfix+'.csv', 'w') as f:
        f.writelines(str(t) + '\n')


def import_results(file):
    d = np.loadtxt(file, dtype='str')
    return d.tolist()


if __name__ == '__main__':
    list = import_results(input_file)

    grouped_dict = group_by(list, column_names[0], postfixes[0])

    for k,v in grouped_dict.iteritems():
        dump_to_files(k, v)
