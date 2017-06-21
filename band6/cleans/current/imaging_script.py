"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess
from glob import glob
import numpy as np

# Create list of measurement sets to use in clean:
# files = glob("../../data_files/*.timing.reweighted.ms")

#Clean variables to be changed
dates = ['../../final_visibilities/aumic_mar_allspws_corrected.ms',
         '../../final_visibilities/aumic_aug_allspws_corrected.ms',
         '../../final_visibilities/aumic_jun_noflare_allspws_corrected.ms']
imsize = 512
pixsize ='0.03arcsec',



visibilities = dates[:]
filename = 'aumic_all'

natural = False
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
all_natural_rms = 1.4494822607957758e-05
all_briggs_rms = 1.578679439262487e-05
all_taper_rms = 1.89688053069e-05
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
        
        
# Uniform
if uniform:
    #Dirty clean to get rms and check if mask is good
    image = filename + "_dirty_uniform"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='uniform',
           niter=0)
    dirty_uniform_rms = imstat(imagename='{}.residual'.format(image),
                region='aumic_rms.region', listit=False)['rms'][0]
    if view:
        #Show dirty image, then clean up and delete all dirty clean files
        viewer(infile=image + '.residual', displaytype='contour')
        print 'uniform dirty rms: {}'.format(dirty_uniform_rms)
        raw_input('mask ready? ')
    subprocess.call("rm -rf {}.*".format(image), shell=True)

    # Clean with correct mask and rms
    image = filename + "_uniform"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='uniform',
           niter=1000000000,
           threshold=dirty_uniform_rms / 2.,
           usemask='user',
           mask=uniform_mask,
           pbmask=None)

    if view: viewer(infile=image + '.image')
    uniform_rms = imstat(imagename='{}.image'.format(image), region='aumic_rms.region',
        listit=False)['rms'][0]
    print 'Uniform clean rms is', uniform_rms

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
        fitsimage='{}.fits'.format(image))
        
# briggs
if briggs:
    #Dirty clean to get rms and check if mask is good
    image = filename + "_dirty_briggs"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='uniform',
           niter=0)
    dirty_briggs_rms = imstat(imagename='{}.residual'.format(image),
                region='aumic_rms.region', listit=False)['rms'][0]
    if view:
        #Show dirty image, then clean up and delete all dirty clean files
        print 'dirty briggs: {}'.format(dirty_briggs_rms)
        viewer(infile=image + '.residual', displaytype='contour')
        raw_input('mask ready? ')
    subprocess.call("rm -rf {}.*".format(image), shell=True)

    # Clean with correct mask and rms
    image = filename + "_briggs"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=visibilities,
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='briggs',
           robust=0.5,
           niter=1000000000,
           threshold=dirty_briggs_rms / 2.,
           usemask='user',
           mask=briggs_mask,
           pbmask=None)

    briggs_rms = imstat(imagename='{}.image'.format(image), region='aumic_rms.region',
        listit=False)['rms'][0]
    print 'briggs clean rms is', briggs_rms
    if view: viewer(infile=image + '.image')

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
        fitsimage='{}.fits'.format(image))
