"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess
from glob import glob
import numpy as np

# Create list of measurement sets to use in clean:
# files = glob("../../data_files/*.timing.reweighted.ms")

#Clean variables to be changed
imsize = 512
pixsize ='0.03arcsec',



visibilities = ['../calibrated_final.ms']
filename = 'aumic_band9'

natural = True
taper = False

natural_mask = 'aumic_band9_mask_natural.region'
taper_mask = 'aumic_band9_mask_taper.region'

view=False

# Concat before clean to fix proper motion offset;
#pointing center of first chronological date is used
if len(visibilities) > 1:
    print('Files to be concatenated are:')
    print(np.array(visibilities))
    print('They will be cleaned under name "{}"'.format(filename))

    raw_input('ok?: ')
    print('Concatenating...')
    subprocess.call("rm -rf {}".format(filename + ".ms"), shell=True)
    concat(vis=visibilities, concatvis=(filename + ".ms"), dirtol='2arcsec' )
    visibilities = filename + ".ms"
else:
    print('Files to be cleaned are {}, under name {}.blah'.format(visibilities, filename))
    raw_input('ok?: ')



#==========================================================
#rms values for centered cleans
natural_rms = 1.2078784e-04
taper_rms =   1.8785533e-04
#==========================================================


print('Cleaning...')
# Natural no taper: residual is the dirty image
if natural:

    #Dirty clean to get rms and check if mask is good
    image = filename+"_dirty_natural"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           niter=0)
    dirty_natural_rms = imstat(imagename='{}.image'.format(image),
        region='aumic_band9_rms.region', listit=False)['rms'][0]

    if view:
        #Show dirty image, then clean up and delete all dirty clean files
        print 'dirty natural rms: {}'.format(dirty_natural_rms)
        viewer(infile=image + '.image', displaytype='contour')
        raw_input('mask ready? ')
    subprocess.call("rm -rf {}.*".format(image), shell=True)

    # Clean with correct mask and rms
    image = filename + "_natural"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           niter=100000000,
           threshold=dirty_natural_rms / 2.,
           usemask='user',
           mask=natural_mask,
           pbmask=None)

    natural_rms = imstat(imagename='{}.image'.format(image),
        region='aumic_band9_rms.region', listit=False)['rms'][0]
    print 'Natural clean rms is', natural_rms
    if view: viewer(infile=image + '.image')

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
               fitsimage='{}.fits'.format(image))





# Natural 200klambda taper:
if taper:
    #Dirty clean to get rms and check if mask is good
    image = filename + "_dirty_taper"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           uvtaper=['200klambda'],
           niter=0)
    dirty_taper_rms = imstat(imagename='{}.image'.format(image),
                region='aumic_band9_rms.region', listit=False)['rms'][0]
    if view:
        #Show dirty image, then clean up and delete all dirty clean files
        print 'dirty taper rms: {}'.format(dirty_taper_rms)
        viewer(infile=image + '.image', displaytype='contour')
        raw_input('mask ready? ')
    subprocess.call("rm -rf {}.*".format(image), shell=True)

    # Clean with correct mask and rms
    image = filename + "_taper"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           uvtaper=['200klambda'],
           niter=1000000000,
           threshold=dirty_taper_rms / 2.,
           usemask='user',
           mask=taper_mask,
           pbmask=None)

    if view: viewer(infile=image + '.image')
    taper_rms = imstat(imagename='{}.image'.format(image), region='aumic_band9_rms.region',
        listit=False)['rms'][0]
    print 'Taper clean rms is', taper_rms

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
        fitsimage='{}.fits'.format(image))
