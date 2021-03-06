from __future__ import division
import numpy as np
import time
import os

from numpy.fft import *


N = 100000
ITERS = 10

if os.getenv('N_ELEMS'):
    N = int(os.getenv('N_ELEMS'))

if os.getenv('N_ITERS'):
    ITERS = int(os.getenv('N_ITERS'))

print "Number of Iters: %d" % ITERS
print "Number of Elems/Iteration: %d" % N
print "\n"

v = np.zeros(N) + 0j
sum = 0.0
elapsed = 0.0

for iter in range(ITERS):
    
    for i in range(len(v)):
        v[i] = (i * 1.0) / N + ((i * 1.0) / N) * 1j

    start = time.time()

    out = fft(v)

    end = time.time()
    elapsed += end - start

    sum += np.sum(np.abs(out))

print "FFT of %d elems" % len(v)
print "Elapsed Time : %.4f" % elapsed, " s"
print "Throughput : %.3f" % ((N*ITERS*16)/(elapsed *1000000)), "MB/s"
print "Sum : %.2e" % (sum/ITERS)
print "\n"
