from python.twostream import two_stream
import numpy as np
import time

tau, w0 = np.loadtxt('tau_and_w0.txt').T
u0 = 0.6427876096865394
Rsfc = 0.25

amean = two_stream(tau, w0, u0, Rsfc)

# write to file
np.savetxt('results/amean_py.txt',amean)

# timing

nt = 10000
start = time.time()
for i in range(nt):
    amean = two_stream(tau, w0, u0, Rsfc)
end = time.time()
print('One run in Python takes',(end-start)/nt,'seconds')