import numpy as np
import matplotlib.pyplot as plt
import sys

data = np.loadtxt(sys.argv[1])
# x, y = np.histogram(data, bins=100)
# print x.shape
# print y.shape
plt.hist(data, bins=100)
plt.show()
