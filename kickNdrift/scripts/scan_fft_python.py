import subprocess
import os
import time
#import signal

print "\nPython simulation\n"

project_dir = '/afs/cern.ch/work/k/kiliakis/git/cpp-benchmark/rfft/'
exe_dir = project_dir + 'python/'
exe_list = ['rfft-bench', 'irfft-bench', 'fft-bench', 'ifft-bench']
#datafiles = '/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-python/__TEST_CASES/main_files/'
outfiles = '/afs/cern.ch/work/k/kiliakis/results/raw/python/fft/v1/'

os.chdir(exe_dir)


n_iterations_list = ['100']
#['100000', '500000', '1000000']
n_points_list = ['1000', '10000', '100000', '500000', '1000000']
#n_points_list = ['100000']
#['10000', '20000', '50000', '100000']
n_threads_list = ['1'] #['1', '2', '4', '8', '10', '14', '18', '22', '28', '32', '56']
repeats = 5
os.chdir(exe_dir)
total_sims = len(n_iterations_list) * len(n_points_list) * \
    len(n_threads_list) * len(exe_list) * repeats


current_sim = 0
for n_iterations in n_iterations_list:
    os.environ['N_ITERS'] = n_iterations
    for n_points in n_points_list:
        os.environ['N_ELEMS'] = n_points
        for n_threads in n_threads_list:
            for exe in exe_list:
                name = exe+'n_p' + n_points + 'n_i' + n_iterations +\
                    'n_thr' + n_threads
                if not os.path.exists(outfiles):
                    os.makedirs(outfiles)
                stdout = open(outfiles + name+'.stdout', 'w')
                stderr = open(outfiles + name+'.stderr', 'w')
                for i in range(0, repeats):
                    print n_iterations, n_points, n_threads, exe, i
                    #start = time.time()
                    subprocess.call(['python', exe+'.py'],
                                    stdout=stdout,
                                    stderr=stdout,
                                    env=os.environ.copy()
                                    )
                    #end = time.time()
                    current_sim += 1
                    #res.write(str(end-start)+'\n')
                    print "%lf %% is completed" % (100.0 * current_sim / total_sims)

