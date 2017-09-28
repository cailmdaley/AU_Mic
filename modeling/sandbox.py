from astrocail import mcmc, plotting
from run6_variable_starflux import make_best_fits
import aumic_fitting
import subprocess as sp
%reload_ext autoreload
%autoreload 2


# run = mcmc.MCMCrun('run6', path='~/Documents/Research/AU_Mic/modeling/run6/run6_chain', nwalkers=50, burn_in=1500); 
run = mcmc.MCMCrun('run6', nwalkers=50, burn_in=1500); aumic_fitting.label_fix(run)

make_best_fits(run)
