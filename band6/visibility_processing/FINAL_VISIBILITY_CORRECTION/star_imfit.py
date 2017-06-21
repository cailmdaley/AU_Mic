
"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess
from glob import glob

# Create list of measurement sets to use in clean:
# files = glob("../../data_files/*.timing.reweighted.ms")

#Clean variables to be changed
dates = ['aumic_aug_centered.concat.ms',
         'aumic_jun_timing_centered.concat.ms',
         'aumic_mar_centered.concat.ms']
imsize = 512
pixsize ='0.03arcsec',


visname = 'aumic_mar.concat.ms'
filename = 'marl_allspws.concat'

natural = True
imfit=False

view=True

print('Cleaning...')
# Natural no taper: residual is the dirty image
if natural:
    image = filename+"_dirty_natural"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visname,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           niter=0)
    dirty_natural_rms = imstat(imagename='{}.residual'.format(image),
                 region='aumic_rms.region', listit=False)['rms'][0]
    # User mask:
    image = filename + "_natural"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visname,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           niter=100000,
           threshold=dirty_natural_rms / 2.,
           usemask='user',
           mask='',
           pbmask=None)
    if view:
        viewer(infile=image + '.image')
    natural_rms = imstat(imagename='{}.image'.format(image),
                 region='aumic_rms.region', listit=False)['rms'][0]
    print 'Natural clean rms is', natural_rms

#Call imfit on region around star
if imfit:
    region = input('Please enter imfit region filename: ')
    fit = imfit(imagename = image+'.image', region = region)
