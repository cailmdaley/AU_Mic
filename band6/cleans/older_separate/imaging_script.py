"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess

# Create list of measurement sets to use in clean:
from glob import glob
files = glob("../data_files/*.ms")
files = np.append(files[4:8], files[9:12])
print(files)
mask='../aumic_larger.mask'

# Rest wavelength = 1351 microns
# ALMA Atenna diameter = 12 m
# Primary beam = lambda/D --> ~23"
# AU Mic ~ 8" across-- going with mask 1/2 of pb


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Concat before clean to try to fix offset caused by proper motion


# concatvis = "aumic_marjune_concat.ms"
# subprocess.call("rm -rf {}".format(concatvis), shell=True)
# concat(vis=files, concatvis=concatvis, dirtol='2arcsec')

# #Natural no taper: residual is the dirty image
# image = "aumic_marjune_concat_dirty_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=concatvis,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# print(rms)
# # User mask:
# image = "aumic_marjune_concat_usermask_natural"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=concatvis,
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
# print(rms) #rms = 1.77061756403e-05
# exportfits(imagename='{}.image'.format(image), fitsimage='{}.fits'.format(image))

# Natural 200klambda taper: residual is the dirty image
image = "aumic_marjune_concat_dirty_200klam"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=concatvis,
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       uvtaper=['200klambda'],
       niter=0)
rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# User mask:
image = "aumic_marjune_concat_usermask_200klam"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=concatvis,
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
listit=False)['rms'][0] #rms= 2.19239682338e-05
print(rms)
exportfits(imagename='{}.image'.format(image), fitsimage='{}.fits'.format(image))
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


# Natural no taper: residual is the dirty image
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
# print(rms) #rms = 1.79587623279e-05
# exportfits(imagename='{}.image'.format(image), fitsimage='{}.fits'.format(image))

# # Natural 200klambda taper: residual is the dirty image
# image = "aumic_marjune_dirty_200klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=files,
#        imagename=image,
#        imsize=512,
#        cell='0.03arcsec',
#        weighting='natural',
#        uvtaper=['200klambda'],
#        niter=0)
# rms = imstat(imagename='{}.residual'.format(image), region='../rms.region', listit=False)['rms'][0]
# # User mask:
# image = "aumic_marjune_usermask_200klam"
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
# rms = imstat(imagename='{}.image'.format(image), region='../rms.region',
# listit=False)['rms'][0] #rms= 2.2775318939238787e-05
# exportfits(imagename='{}.image'.format(image), fitsimage='{}.fits'.format(image))
