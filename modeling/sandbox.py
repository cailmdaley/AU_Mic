from astrocail import mcmc, plotting
import aumic_fitting
import subprocess as sp
%reload_ext autoreload
%autoreload 2


run = mcmc.MCMCrun('run6', path='~/Documents/Research/AU_Mic/modeling/run6/run6_chain', nwalkers=50, burn_in=1500); 

aumic_fitting.label_fix(run)

run.converged


# run.kde()
