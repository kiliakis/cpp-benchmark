from __future__ import division
import numpy as np
import time
import os


N = 100000
ITERS = 1

if os.getenv('N_ELEMS'):
    N = int(os.getenv('N_ELEMS'))

if os.getenv('N_ITERS'):
    ITERS = int(os.getenv('N_ITERS'))

print "Number of turns: %d" % ITERS
print "Number of points: %d" % N
print "\n"

a = np.random.rand(N)
b = np.random.rand(N)
sum = 0.0
elapsed = 0.0

start = time.time()
for iter in range(ITERS):

    c = np.minimum(a, b)

elapsed = time.time() - start
sum = np.sum(c)

throughput = N * ITERS / elapsed / 1.e6
print "c = min(a, b) bench"
print "Elapsed Time : %.4f sec" % elapsed
print "Throughput : %.3f MB/sec" % throughput
print "Sum : %.5e" % sum
print "\n"
