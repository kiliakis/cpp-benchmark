import subprocess
import os
import time
#import signal


print('\nCpp simulation\n')

# os.environ['GOMP_CPU_AFFINITY'] = '0 28 1 29 2 30 3 31 4 32 5 33 6 34 7 35 8 36 9 37 10 ' + \
#    '38 11 39 12 40 13 41 14 42 15 43 16 44 17 45 18 46 19 47 20 48 21 49 22 50 23 51 24 52 25 53 26 54 27 55'
this_directory = os.path.dirname(os.path.realpath(__file__)) + "/"

project_dir = this_directory + '../'
exe_dir = project_dir + 'build/exe/'
exe_list = [
    'smooth_histo1',
    'smooth_histo2',
    'smooth_histo3',
    'smooth_histo4',
    ]


#datafiles = '/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/Release/'
outfiles = project_dir + 'results/raw/perftest/'

n_iterations_list = ['1000']
n_points_list = ['1000000', '2000000', '4000000', '8000000', '16000000']
# n__list = ['1000000', '2000000', '4000000', '8000000', '16000000']
# n_points_list = ['1000000']
# n_points_list = ['16000000']
n_threads_list = ['1', '4', '8', '14', '28']
# n_threads_list = ['4']

repeats = 4
os.chdir(exe_dir)
total_sims = len(n_iterations_list) * len(n_points_list) * \
    len(n_threads_list) * len(exe_list) * repeats

current_sim = 0
for n_iterations in n_iterations_list:
    for n_points in n_points_list:
        n_slices = str(int(n_points) // 1000)
        for n_threads in n_threads_list: 
            os.environ['OMP_NUM_THREADS'] = str(n_threads)
            for exe in exe_list:
                name = '{}-n_p{}-n_s{}-n_i{}-n_thr{}'.format(exe, n_points, n_slices, n_iterations, n_threads)
                if not os.path.exists(outfiles):
                    os.makedirs(outfiles)
                #res = open(outfiles + name+'.res', 'w')
                stdout = open(outfiles + name+'.stdout', 'w')
                stderr = open(outfiles + name+'.stderr', 'w')
                for i in range(0, repeats):
                    exe_args = ['./'+exe,
                                '-p' + n_points,
                                '-m' + n_threads,
                                '-t' + n_iterations,
                                '-s' + n_slices
                                ]
                    print(name, i)
                    #start = time.time()
                    subprocess.call(exe_args,
                                    stdout=stdout,
                                    stderr=stderr,
                                    env=os.environ.copy()
                                    )
                    #end = time.time()
                    current_sim += 1
                    # res.write(str(end-start)+'\n')
                    print("%.2f %% is completed" % (100.0 * current_sim / total_sims))
