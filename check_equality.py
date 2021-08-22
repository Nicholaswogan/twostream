import numpy as np
amean_cpp = np.loadtxt('results/amean_cpp.txt').T
amean_f90 = np.loadtxt('results/amean_fortran.txt').T
amean_py = np.loadtxt('results/amean_py.txt').T
print(np.all(np.isclose(amean_cpp,amean_f90)))
print(np.all(np.isclose(amean_cpp,amean_py)))