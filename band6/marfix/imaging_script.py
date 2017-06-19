"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess

imsize = 512
pixsize ='0.03arcsec',



visibilities = ['aumic_mar_allspws.concat.ms']
filename = 'aumic_mar_with_star'
natural = False
natural_mask = 'aumic_jun_mask.region'
taper = True
taper_mask = 'aumic_mar_mask.region'
view=True

# Concat before clean to fix proper motion offset;
#pointing center of first chronological date is used
if len(visibilities) > 1:
    print('Files to be cleaned are {}, under name {}.concat.ms'.format(visibilities, filename))
    raw_input('ok?: ')
    print('Concatenating...')
    subprocess.call("rm -rf {}".format(filename + ".concat.ms"), shell=True)
    concat(vis=visibilities, concatvis=(filename + ".concat.ms"), dirtol='2arcsec' )
else:
    print('Files to be cleaned are {}, under name {}.blah'.format(visibilities, filename))
    raw_input('ok?: ')

#==========================================================
#rms values 
rms = 2.96353882732e-05
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
        viewer(infile=image + '.residual', displaytype='contour')
        print 'dirty natural: {}'.format(dirty_natural_rms)
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

    if view: viewer(infile=image + '.image')
    natural_rms = imstat(imagename='{}.image'.format(image),
        region='aumic_rms.region', listit=False)['rms'][0]
    print 'Natural clean rms is', natural_rms

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
        viewer(infile=image + '.residual', displaytype='contour')
        print 'dirty taper: {}'.format(dirty_taper_rms)
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
