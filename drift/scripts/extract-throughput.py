#!/usr/bin/python
import os
import csv
import statistics as stat
import sys
import numpy as np
import re
# from prettytable import PrettyTable

       
def string_between(string, before, after):
    temp = string
    if before:
        temp = temp.split(before)[1]
    if after:
        temp = temp.split(after)[0]
    return temp


def extract_results(input, out):
    header = ['bench', 'n_threads', 'n_points', 'n_iterations', 'alpha',
              'time(s)', 'throughput(mp/sec)', 'throughput/th', 't_std', 'th_std', 'th/th_std']
    records = []
    regexp = 'Number of turns:\s(\d+)\sNumber of points:\s(\d+)\sNumber of openmp threads:\s(\d+)\sAlpha order:\s(\d+)\sDrift\s(\w+)\sElapsed time:\s(.*)\ssec\sThroughput:\s(.*)\sMP/sec\sThroughput/thread:\s(.*)\sMP/sec'
    regexp = re.compile(regexp)
    for dirs, subdirs, files in os.walk(input):
        for file in files:
            if ('.stdout' not in file):
                continue
            # alpha = string_between(file, 'n_a', 'n_i')
            # points = string_between(file, 'n_p', 'n_i')
            # iters = string_between(file, 'n_i', 'n_thr')
            # exe = string_between(file, '', '_np')
            # #implementation = dirs.split('/')[-2]
            lines = open(os.path.join(dirs, file), 'r').read()
            match = regexp.findall(lines)
            # print(len(match), match[0])
            if match:
                assert(len(match[0]) == 8)
                l = []
                turns = match[0][0]
                points = match[0][1]
                threads = match[0][2]
                alpha = match[0][3]
                bench = match[0][4]
                for m in match:
                    l.append(m[5:])
                l = np.array(l, float)
                avg = np.round(np.mean(l, axis=0), 3)
                std = np.round(np.std(l, axis=0), 3)
                records.append([bench, threads, points, turns, alpha] + list(avg) + list(std))

            # for line in open(os.path.join(dirs, file), 'r'):
            #     if 'Elapsed' in line:
            #         time = string_between()
            #     if 'Throughput/thread' in line:

            #     if 'Throughput' in line:
            #         temp = string_between(line, ':', 'M')
            #         l.append(float(temp.strip()))
            # if l:
            #     records.append([exe, threads, points, iters,
            #                     stat.mean(l), stat.stdev(l)])
            #     print(file)
            #     percent = 100.0 * stat.stdev(l) / stat.mean(l)
            #     if percent > 10:
            #         print("The previous file has %.2f %% error" % percent)
            # else:
            #     print("I found an empty file called %s" % file)
    records.sort(key=lambda a: (a[0], int(a[1]), int(a[2]), int(a[3]), int(a[4])))
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
    d = np.loadtxt(output_file, dtype='str')
    return d.tolist()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("You should specify input directory and output file")
        exit(-1)
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    extract_results(input_dir, output_file)
    
    
