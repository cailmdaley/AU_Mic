"""
Contains various cleaning methods for AU Mic--to be run in aumic_composite_cleans directory.
"""
import subprocess

# Create list of measurement sets to use in clean:
from glob import glob
files = glob("../data_files/*.ms")
mask='../aumic.mask'

# Rest wavelength = 1351 microns
# ALMA Atenna diameter = 12 m
# Primary beam = lambda/D --> ~23"
# AU Mic ~ 8" across-- going with mask 1/2 of pb
# exportfits(imagename=image+'.image', fitsimage=image + '.fits')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auto-thresh clean
# image = "aumic_auto-thresh_briggs2_niter10000"
# subprocess.call("rm -rf {}*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=2,
#        niter=10000,
#        usemask='auto-thresh',
#        pbmask=0.5)
#
# viewer(infile=image + '.image')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Natural no taper: residual is the dirty image
image = "aumic_dirty_natural"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=files,
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       niter=0)
rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# User mask:
image = "aumic_usermask_natural"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=files,
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       niter=100000,
       threshold=rms/2.,
       usemask='user',
       mask = mask,
       pbmask=None)
viewer(infile=image+'.image')
rms = imstat(imagename='{}.image'.format(image), region='rms.region', listit=False)['rms'][0] # rms=1.4716037185280584e-05


# # Natural 100klambda taper: residual is the dirty image
# image = "aumic_dirty_natural_100klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['100klambda'],
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_usermask_natural_100klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['100klambda'],
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')


# # Natural 200klambda taper: residual is the dirty image
# image = "aumic_dirty_natural_200klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['200klambda'],
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_usermask_natural_200klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['200klambda'],
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')
# rms = imstat(imagename='{}.image'.format(image), region='rms.region',
# listit=False)['rms'][0] #rms=1.9399181837798096e-05

# # Natural 300klambda taper: residual is the dirty image
# image = "aumic_dirty_natural_300klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['300klambda'],
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_usermask_natural_300klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['300klambda'],
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')


# # Natural 400klambda taper: residual is the dirty image
# image = "aumic_dirty_natural_400klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['400klambda'],
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_usermask_natural_400klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['400klambda'],
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')


# # Natural 500klambda taper: residual is the dirty image
# image = "aumic_dirty_natural_500klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['500klambda'],
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_usermask_natural_500klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['500klambda'],
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Briggs Robust = 2: residual is the dirty image
# image = "aumic_dirty_briggs2"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=2,
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_usermask_briggs2"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=2,
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Briggs Robust = 1: residual is the dirty image
# image = "aumic_dirty_briggs1"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=1,
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]

# User mask:
# image = "aumic_usermask_briggs1"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=1,
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Briggs Robust = 0.5: residual is the dirty image
# image = "aumic_dirty_briggs0.5"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=0.5,
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
#
# # User mask:
# image = "aumic_usermask_briggs0.5"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=0.5,
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Briggs Robust = 0: residual is the dirty image
# image = "aumic_dirty_briggs0"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=0,
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]

# User mask:
# image = "aumic_usermask_briggs0"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=0,
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #Briggs Robust = -2: residual is the dirty image
# image = "aumic_dirty_briggs-2"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=-2,
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='rms.region', listit=False)['rms'][0]
#
# # User mask:
# image = "aumic_usermask_briggs-2"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='briggs',
#        robust=-2,
#        niter=100000,
#        threshold=rms/2.,
#        usemask='user',
#        mask = mask,
#        pbmask=None)
# viewer(infile=image + '.image')
