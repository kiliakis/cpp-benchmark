synch_rad1
std, non parallel dist gen, separate loops

synch_rad3
std, parallel dist gen, separate loops

synch_rad4
std, parallel dist gen, single loop

synch_rad7
std, parallel dist gen, single loop, vectorized

synch_rad2
boost, non parallel dist gen, separate loops

synch_rad5
boost, parallel dist gen, separate loops

synch_rad6
boost, parallel dist gen, single loop

synch_rad8
boost, parallel dist gen, single loop, vectorized

synch_rad9
boost, parallel dist gen, single loop, vectorized, only one mem alloc per thread


Which are interesting to run?
all of them 
Sizes? N_t 100, runs 5, n_p 1, 2, 4, 8, 16, 32, 64
threads? 1, 4, 8, 14, 28

Plots
What plots do I need?
at first a simple table to compare, single thread, particles 16M:
all implementations

Then a plot to compare scalability:
of all boost implementations, and all std implementations

Okay I can use the one with 16M particles + the table