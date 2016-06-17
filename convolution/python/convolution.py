from __future__ import division
import numpy as np
import time
import os


N = 100000
ITERS = 10

if os.getenv('N_ELEMS'):
    N = int(os.getenv('N_ELEMS'))

if os.getenv('N_ITERS'):
    ITERS = int(os.getenv('N_ITERS'))

print "Number of Iters: %d" % ITERS
print "Number of Elems/Iteration: %d" % N
print "\n"

v = np.zeros(N)
kernel = np.zeros(N/2)
sum = 0.0
elapsed = 0.0

for iter in range(ITERS):
    
    for i in range(len(v)):
        v[i] = (i * 1.0) / N

    for i in range(len(kernel)):
        kernel[i] = (i * 1.0) / (10.0 * N)

    start = time.time()

    out = np.convolve(v, kernel)

    end = time.time()
    elapsed += end - start

    sum += np.sum(np.abs(out))

print "Convolution of %d elems" % len(v)
print "Elapsed Time : %.4f" % elapsed, " s"
print "Throughput : %.3f" % ((N*ITERS*8)/(elapsed *1000000)), "MB/s"
print "Sum : %.2e" % (sum/ITERS)
print "\n"
