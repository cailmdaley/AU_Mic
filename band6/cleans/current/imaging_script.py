"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess
from glob import glob
import numpy as np

# Create list of measurement sets to use in clean:
# files = glob("../../data_files/*.timing.reweighted.ms")

#Clean variables to be changed
dates = ['../../final_visibilities/aumic_mar_allspws_FINAL.ms',
         '../../final_visibilities/aumic_aug_allspws_FINAL.ms',
         '../../final_visibilities/aumic_jun_allspws_FINAL.ms']
dates2 = [
    '../../visibility_processing/FINAL_VISIBILITY_CORRECTION/aumic_mar_allspws.fixvis.uvsub.ms',
    '../../visibility_processing/FINAL_VISIBILITY_CORRECTION/aumic_aug_allspws.fixvis.uvsub.ms',
    '../../visibility_processing/FINAL_VISIBILITY_CORRECTION/aumic_jun_noflare_allspws.fixvis.uvsub.ms']
imsize = 512
pixsize ='0.03arcsec',



visibilities = dates
filename = 'aumic_all_unweighted'

natural = True
taper = True
uniform = False
briggs = False

natural_mask = 'aumic_mask_all_natural.region'
taper_mask = 'aumic_mask_all_taper.region'
uniform_mask = 'aumic_mask_all_uniform.region'
briggs_mask = 'aumic_mask_all_briggs.region'

view=True

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
all_natural_rms = 1.57988852152e-05
all_taper_rms = 1.97743011086e-05
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
    dirty_natural_rms = imstat(imagename='{}.residual'.format(image),
        region='aumic_rms.region', listit=False)['rms'][0]

    if view:
        #Show dirty image, then clean up and delete all dirty clean files
        print 'dirty natural: {}'.format(dirty_natural_rms)
        viewer(infile=image + '.residual', displaytype='contour')
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
        region='aumic_rms.region', listit=False)['rms'][0]
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
    dirty_taper_rms = imstat(imagename='{}.residual'.format(image),
                region='aumic_rms.region', listit=False)['rms'][0]
    if view:
        #Show dirty image, then clean up and delete all dirty clean files
        print 'dirty taper: {}'.format(dirty_taper_rms)
        viewer(infile=image + '.residual', displaytype='contour')
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
    taper_rms = imstat(imagename='{}.image'.format(image), region='aumic_rms.region',
        listit=False)['rms'][0]
    print 'Taper clean rms is', taper_rms

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
        fitsimage='{}.fits'.format(image))
        
