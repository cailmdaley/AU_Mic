import pyximport; pyximport.install()

from disk_model import test, test_py


cydisk = test.Disk()
pydisk = test_py.Disk()


import cProfile
cProfile.run('pydisk.set_structure()')

%timeit pydisk.set_structure()
10 loops, best of 3: 30.9 ms per loop

%timeit cydisk.set_structure() # copy-and-pasted python
10 loops, best of 3: 32.6 ms per loop

%timeit cydisk.set_structure() # cython loop
10 loops, best of 3: 23.2 ms per loop

%timeit cydisk.set_structure() # cdef class
10 loops, best of 3: 31.5 ms per loop

%timeit cydisk.set_structure() # cdef class and cpdef func
10 loops, best of 3: 29.8 ms per loop











# from astrocail import mcmc, plotting
# from run6_variable_starflux import make_best_fits
# import aumic_fitting
# import subprocess as sp
# %reload_ext autoreload
# %autoreload 2
# 
# 
# # run = mcmc.MCMCrun('run6', path='~/Documents/Research/AU_Mic/modeling/run6/run6_chain', nwalkers=50, burn_in=1500); 
# run = mcmc.MCMCrun('run6', nwalkers=50, burn_in=1500); aumic_fitting.label_fix(run)
# 
# make_best_fits(run)
