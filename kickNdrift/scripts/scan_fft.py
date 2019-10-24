import subprocess
import os
import time
#import signal


print '\nCpp simulation\n'

#os.environ['GOMP_CPU_AFFINITY'] = '0 28 1 29 2 30 3 31 4 32 5 33 6 34 7 35 8 36 9 37 10 ' + \
#    '38 11 39 12 40 13 41 14 42 15 43 16 44 17 45 18 46 19 47 20 48 21 49 22 50 23 51 24 52 25 53 26 54 27 55'

project_dir = '/afs/cern.ch/work/k/kiliakis/git/cpp-benchmark/rfft/'
exe_dir = project_dir + 'build-2/'
exe_list = ['./irfft-bench', './rfft-bench']#, './fft-bench', './ifft-bench']
#datafiles = '/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/Release/'
outfiles = '/afs/cern.ch/work/k/kiliakis/results/raw/cpp/fft/v8/'

n_iterations_list = ['100']
#['100000', '500000', '1000000']
n_points_list = ['1000000']#['1000', '10000', '100000', '500000', '1000000']
#n_points_list = ['100000']
#['10000', '20000', '50000', '100000']
n_threads_list = ['50'] #['1', '2', '4', '8', '10', '14', '18', '22', '28', '32', '56']
repeats = 5
os.chdir(exe_dir)
total_sims = len(n_iterations_list) * len(n_points_list) * \
    len(n_threads_list) * len(exe_list) * repeats

current_sim = 0
for n_iterations in n_iterations_list:
    for n_points in n_points_list:
        for n_threads in n_threads_list:
            for exe in exe_list:
                name = exe+'n_p' + n_points + 'n_i' + n_iterations +\
                    'n_thr' + n_threads
                if not os.path.exists(outfiles):
                    os.makedirs(outfiles)
                #res = open(outfiles + name+'.res', 'w')
                stdout = open(outfiles + name+'.stdout', 'w')
                stderr = open(outfiles + name+'.stderr', 'w')
                for i in range(0, repeats):
                    exe_args = [exe,
                                '-n' + n_points,
                                '-t' + n_threads,
                                '-i' + n_iterations
                                ]
                    print n_iterations, n_points, n_threads, exe, i
                    #start = time.time()
                    subprocess.call(exe_args,
                                    stdout=stdout,
                                    stderr=stderr
                                    )
                    #end = time.time()
                    current_sim += 1
                    #res.write(str(end-start)+'\n')
                    print "%.2f %% is completed" % (100.0 * current_sim / total_sims)
