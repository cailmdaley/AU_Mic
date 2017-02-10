"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess

# Create list of measurement sets to use in clean:
from glob import glob
files = glob("../data_files/*.ms")
files = files[4:]
print(files)
mask='../aumic_larger.mask'

# Rest wavelength = 1351 microns
# ALMA Atenna diameter = 12 m
# Primary beam = lambda/D --> ~23"
# AU Mic ~ 8" across-- going with mask 1/2 of pb
# exportfits(imagename=image+'.image', fitsimage=image + '.fits')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Natural no taper: residual is the dirty image
# image = "aumic_18aug_dirty_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# # User mask:
# image = "aumic_18aug_usermask_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image+'.image')
# rms = imstat(imagename='{}.image'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# #18aug rms=2.5078490580199286e-05
# #24jun rms=2.02323226404e-05aumic_18aug_usermask_natural.


# # Natural no taper: residual is the dirty image
# image = "aumic_26mar_dirty_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# # User mask:
# image = "aumic_26mar_usermask_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image+'.image')
# rms = imstat(imagename='{}.image'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# #rms = 2.60297147179e-05


# # Natural no taper: residual is the dirty image
# image = "aumic_marjune_dirty_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# # User mask:
# image = "aumic_marjune_usermask_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image+'.image')
# rms = imstat(imagename='{}.image'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# #rms = 1.7922860934049822e-05

# Natural 200klambda taper: residual is the dirty image
image = "aumic_marjune_dirty_200klam"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=files,
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       uvtaper=['200klambda'],
       niter=0)
rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# User mask:
image = "aumic_marjune_usermask_200klam"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=files,
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       uvtaper=['200klambda'],
       niter=100000,
       threshold=rms/2.,
       usemask='user',
       mask = mask,
       pbmask=None)
viewer(infile=image + '.image')
rms = imstat(imagename='{}.image'.format(image), region='../rms.region',
listit=False)['rms'][0] #rms=2.2735775928595103e-05
