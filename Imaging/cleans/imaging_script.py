"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess
from glob import glob

# Create list of measurement sets to use in clean:
files = glob("../../data_files/*.ms")

#Clean variables to be changed
concat_files = files[4:8]
filename = 'aumic_jun'
mask = 'aumic_jun_mask.region'
imsize = 512
pixsize ='0.03arcsec',

# Concat before clean to fix proper motion offset
subprocess.call("rm -rf {}".format(filename + ".concat.ms"), shell=True)
concat(vis=concat_files, concatvis=(filename + ".concat.ms"), dirtol='2arcsec' )

#==========================================================
#rms values for various cleans
aug_natural_rms = 3.44968793797e-05
june_natural_rms = 3.37951933034e-05
augjun_natural_rms = 2.9192528018e-05
marjun_natural_rms = 2.27807577176e-05
composite_natural_rms = 2.06453660212e-05
marjun_200klam_rms = 2.99919338431e-05
mar_200klam_rms = 2.92010208796e-05
composite_200klam_rms = 2.7314921681e-05
#==========================================================

# Natural no taper: residual is the dirty image
image = filename+"_dirty_natural"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis= filename + ".concat.ms",
       imagename=image,
       imsize=imsize,
       cell=pixsize,
       weighting='natural',
       niter=0)
rms = imstat(imagename='{}.residual'.format(image),
             region='aumic_rms.region', listit=False)['rms'][0]

# User mask:
image = filename + "_natural"
subprocess.call("rm -rf {}.*".format(image), shell=True)
tclean(vis=filename + ".concat.ms",
       imagename=image,
       imsize=imsize,
       cell=pixsize,
       weighting='natural',
       niter=100000,
       threshold=rms / 2.,
       usemask='user',
       mask=mask,
       pbmask=None)
#viewer(infile=image + '.image')
rms = imstat(imagename='{}.image'.format(image),
             region='aumic_rms.region', listit=False)['rms'][0]
print(rms)

# Export to .fits
exportfits(imagename='{}.image'.format(image),
           fitsimage='{}.fits'.format(image))


# # Natural 200klambda taper: residual is the dirty image
# image = filename + "_dirty_200klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=filename + ".concat.ms",
#       imagename=image,
#       imsize=imsize,
#       cell=pixsize,
#       weighting='natural',
#       uvtaper=['200klambda'],
#       niter=0)
# rms = imstat(imagename='{}.residual'.format(image),
#             region='aumic_rms.region', listit=False)['rms'][0]
#
# # User mask:
# image = filename + "_200klam"
# subprocess.call("rm -rf {}.*".format(image), shell=True)
# tclean(vis=filename + ".concat.ms",
#       imagename=image,
#       imsize=imsize,
#       cell=pixsize,
#       weighting='natural',
#       uvtaper=['200klambda'],
#       niter=100000,
#       threshold=rms / 2.,
#       usemask='user',
#       mask=mask,
#       pbmask=None)
# viewer(infile=image + '.image')
# rms = imstat(imagename='{}.image'.format(image), region='aumic_rms.region',
#             listit=False)['rms'][0]
# print(rms)
#
# # Export to .fits
# exportfits(imagename='{}.image'.format(image),
#           fitsimage='{}.fits'.format(image))
