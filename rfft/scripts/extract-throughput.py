#!/usr/bin/python
import os
import csv
import statistics as stat
import sys
import numpy as np

from prettytable import PrettyTable

       
def string_between(string, before, after):
    temp = string
    if before:
        temp = temp.split(before)[1]
    if after:
        temp = temp.split(after)[0]
    return temp


def extract_results(input, out):
    header = ['fft_type', 'n_threads', 'n_points', 'n_iterations',
              'Throughput', 'Stdev']
    records = []
    for dirs, subdirs, files in os.walk(input):
        for file in files:
            if ('.stderr' not in file):
                continue
            l = []
            threads = string_between(file, 'n_thr', '.')
            points = string_between(file, 'n_p', 'n_i')
            iters = string_between(file, 'n_i', 'n_thr')
            exe = string_between(file, '', '-bench')
            #implementation = dirs.split('/')[-2]
            for line in open(os.path.join(dirs, file), 'r'):
                if line.strip():
                    temp = string_between(line, ':', 'M')
                    l.append(float(temp.strip()))
            if l:
                records.append([exe, threads, points, iters,
                                stat.mean(l), stat.stdev(l)])
                print file
                percent = 100.0 * stat.stdev(l) / stat.mean(l)
                if percent > 10:
                    print "The previous file has %.2f %% error" % percent
            else:
                print "I found an empty file called %s" % file
    records.sort(key=lambda a: (a[0], int(a[1]), int(a[2]), int(a[3])))
    t = PrettyTable(header)
    t.align = 'l'
    t.border = False
    for r in records:
        t.add_row(r)
    with open(out, 'w') as f:
        f.writelines(str(t) + '\n')
    #writer = csv.writer(open(out, 'w'), lineterminator='\n', delimiter='\t')
    #writer.writerow(header)
    #writer.writerows(records)


def import_results(output_file):
    d = np.loadtxt(output_file, dtype='str')
    return d.tolist()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "You should specify input directory and output file"
        exit(-1)
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    extract_results(input_dir, output_file)
    
    
