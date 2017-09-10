from astrocail import mcmc, plotting
import aumic_fitting
import subprocess as sp
%reload_ext autoreload
%autoreload 2


run = mcmc.MCMCrun('run7', path='~/Documents/Research/AU_Mic/modeling/run7/run7_chain', nwalkers=18); 

aumic_fitting.label_fix(run)

run.converged
run.main


# run.kde()
