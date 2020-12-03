#!/usr/bin/python
import os
import csv
import sys
import numpy as np

# from prettytable import PrettyTable

       
def string_between(string, before, after):
    temp = string
    if before:
        temp = temp.split(before)[1]
    if after:
        temp = temp.split(after)[0]
    return temp


def extract_results(input, out):
    header = ['exe', 'n_threads', 'n_points', 'n_slices', 'n_iterations',
              'time', 'stdev']
    records = []
    for dirs, subdirs, files in os.walk(input):
        for file in files:
            if ('.stdout' not in file):
                continue
            l = []
            threads = string_between(file, 'n_thr', '.')
            points = string_between(file, 'n_p', '-')
            iters = string_between(file, 'n_i', '-')
            slices = string_between(file, 'n_s', '-')
            exe = string_between(file, '', '-')
            for line in open(os.path.join(dirs, file), 'r'):
                if 'Elapsed time' in line:
                    temp = string_between(line, ':', 'sec')
                    l.append(float(temp.strip()))
            if l:
                records.append([exe, threads, points, slices, iters,
                                np.mean(l), np.std(l)])
                print(file)
                percent = 100.0 * np.std(l) / np.mean(l)
                if percent > 10:
                    print("The previous file has %.2f %% error" % percent)
            else:
                print("I found an empty file called %s" % file)
    records.sort(key=lambda a: (a[0], int(a[1]), int(a[2]), int(a[3])))
    # t = PrettyTable(header)
    # t.align = 'l'
    # t.border = False
    # for r in records:
    #     t.add_row(r)
    # with open(out, 'w') as f:
    #     f.writelines(str(t) + '\n')
    writer = csv.writer(open(out, 'w'), lineterminator='\n', delimiter='\t')
    writer.writerow(header)
    writer.writerows(records)


def import_results(output_file):
    d = np.genfromtxt(output_file, dtype='str')
    return d.tolist()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("You should specify input directory and output file")
        exit(-1)
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    extract_results(input_dir, output_file)
    
    
