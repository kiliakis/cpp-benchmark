import time
begin = time.time()
start = time.time()
import numpy as np
print "imports: ", time.time() - start
N = 1000000
ITERS = 1

# if os.getenv('N_ELEMS'):
#     N = int(os.getenv('N_ELEMS'))

# if os.getenv('N_ITERS'):
#     ITERS = int(os.getenv('N_ITERS'))

print "Number of turns: %d" % ITERS
print "Number of points: %d" % N
print "\n"

start = time.time()
a = np.random.rand(N)
b = np.random.rand(N)
end = time.time()
print "initialization: ", end - start

sum = 0.0
start = time.time()
for iter in range(ITERS):
    a += b
end = time.time()
elapsed = end - start
print "run time: ", end - start

start = time.time()
sum = np.sum(a)
end = time.time()
print "finalization: ", end - start

throughput = 1. * N * ITERS / elapsed / 1e6
print "a += b bench"
print "Elapsed Time : %.4f sec" % elapsed
print "Throughput : %.3f MB/sec" % throughput
print "Sum : %.5e" % sum
print "\n"
print "Total time: ", time.time() - begin