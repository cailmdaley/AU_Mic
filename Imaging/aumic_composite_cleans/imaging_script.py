"""
Contains various cleaning methods for AU Mic--to be run in aumic_composite_cleans directory.
"""
import subprocess

#Create list of measurement sets to use in clean:
from glob import glob
files = glob("../data_files/*.ms")

#Rest wavelength = 1351 microns
#ALMA Atenna diameter = 12 m
#Primary beam = lambda/D --> ~23"
#AU Mic ~ 8" across-- going with mask 1/2 of pb

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Auto-thresh clean
# image = "aumic_auto-thresh_niter10000"
# subprocess.call("rm -rf {}*".format(image), shell=True)
# tclean(vis=files,
# imagename=image,
# imsize=512,
# cell='0.03arcsec',
# weighting = 'briggs',
# robust = 2,
# niter=10000,
# usemask='auto-thresh',
# pbmask = 0.5)
#
# viewer(infile = image+'.image')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Create dirty map with uniform weighting to get noise
image = "aumic_dirty"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=files,
imagename=image,
imsize=512,
cell='0.03arcsec',
weighting = 'uniform',
niter=1,
usemask='pb',
pbmask = 0.5)

viewer(infile = image+'.image')
