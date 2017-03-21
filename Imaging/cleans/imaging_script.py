"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess

# Create list of measurement sets to use in clean:
from glob import glob
files = glob("../data_files/*.ms")
mask = '../aumic_larger.mask'

# Concat before clean to fix proper motion offset
files = files
concat_name = "aumic_composite"
subprocess.call("rm -rf {}".format(concat_name + "_concat.ms"), shell=True)
concat(vis=files, concatvis=(concat_name + "_concat.ms"), dirtol='2arcsec' )

#==========================================================
#rms values for various cleans
marjune_rms = 0
#==========================================================

# Natural no taper: residual is the dirty image
image = concat_name+"_dirty_natural"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=concat_name + "_concat.ms",
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       niter=0)
rms = imstat(imagename='{}.residual'.format(image),
             region='../rms.region', listit=False)['rms'][0]

# User mask:
image = concat_name + "_usermask_natural"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=concat_name + "_concat.ms",
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       niter=100000,
       threshold=rms / 2.,
       usemask='user',
       mask=mask,
       pbmask=None)
viewer(infile=image + '.image')
rms = imstat(imagename='{}.image'.format(image),
             region='../rms.region', listit=False)['rms'][0]
print(rms)

# Export to .fits
exportfits(imagename='{}.image'.format(image),
           fitsimage='{}.fits'.format(image))


# Natural 200klambda taper: residual is the dirty image
image = concat_name + "_dirty_200klam"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=concat_name + "_concat.ms",
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       uvtaper=['200klambda'],
       niter=0)
rms = imstat(imagename='{}.residual'.format(image),
             region='../rms.region', listit=False)['rms'][0]

# User mask:
image = concat_name + "_usermask_200klam"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=concat_name + "_concat.ms",
       imagename=image,
       imsize=512,
       cell='0.03arcsec',
       weighting='natural',
       uvtaper=['200klambda'],
       niter=100000,
       threshold=rms / 2.,
       usemask='user',
       mask=mask,
       pbmask=None)
viewer(infile=image + '.image')
rms = imstat(imagename='{}.image'.format(image), region='../rms.region',
             listit=False)['rms'][0]
print(rms)

# Export to .fits
exportfits(imagename='{}.image'.format(image),
           fitsimage='{}.fits'.format(image))
