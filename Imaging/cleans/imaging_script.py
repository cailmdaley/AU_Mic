"""
Contains various cleaning methods for AU Mic--to be run in seperate_cleans directory.
"""
import subprocess
from glob import glob

# Create list of measurement sets to use in clean:
# files = glob("../../data_files/*.timing.reweighted.ms")

#Clean variables to be changed
dates = ['aumic_aug.concat.ms',
         'aumic_jun_timing.concat.ms',
         'aumic_mar.concat.ms']
concat_files = ['../../data_files/24jun2015_aumic1_spw1.timing.reweighted.ms']
print 'Files to be cleaned are:', concat_files
raw_input('ok?: ')
# concat_files.append('../Visibility_processing/24jun2015_aumic1_spw3.timeflag.reweighted.ms')

filename = 'aumic_jun_timing_spw1'
imsize = 512
pixsize ='0.03arcsec',


natural = True
natural_mask = 'aumic_jun_mask.region'

taper = False
taper_mask = 'aumic_marjun_short_bsl_mask.region'

view=True

# Concat before clean to fix proper motion offset;
#pointing center of first chronological date is used
# if len(concat_files) > 1:
print('Concatenating...')
subprocess.call("rm -rf {}".format(filename + ".concat.ms"), shell=True)
concat(vis=concat_files, concatvis=(filename + ".concat.ms"), dirtol='2arcsec' )

#==========================================================
#rms values for various cleans
aug_natural_rms = 3.83226033591e-05
jun_timing_rms = 2.79580763163e-05
jun_timing_spw1_rms = 4.26983060606e-05
mar_200klam_rms = 2.92891236313e-05
augjun_natural_rms = 2.25834228331e-05
composite_natural_rms = 1.79467988346e-05
composite_200klam_rms = 2.14701703953e-05
marjun_natural_rms = 2.02876035473e-05
marjun_200klam_rms = 2.39355267695e-05
#==========================================================


print('Cleaning...')
# Natural no taper: residual is the dirty image
if natural:
    image = filename+"_dirty_natural"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis= filename + ".concat.ms",
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           niter=0)
    dirty_natural_rms = imstat(imagename='{}.residual'.format(image),
                 region='aumic_rms.region', listit=False)['rms'][0]
    if view:
        viewer(infile=image + '.residual')
        print 'dirty natural: {}'.format(dirty_natural_rms)

    # User mask:
    image = filename + "_natural"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=filename + ".concat.ms",
           imagename=image,
           imsize=imsize,
           cell=pixsize,
           weighting='natural',
           niter=100000,
           threshold=dirty_natural_rms / 2.,
           usemask='user',
           mask=natural_mask,
           pbmask=None)
    if view:
        viewer(infile=image + '.image')
    natural_rms = imstat(imagename='{}.image'.format(image),
                 region='aumic_rms.region', listit=False)['rms'][0]
    print 'Natural clean rms is', natural_rms

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
               fitsimage='{}.fits'.format(image))

# Natural 200klambda taper: residual is the dirty image
if taper:
    image = filename + "_dirty_200klam"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=filename + ".concat.ms",
          imagename=image,
          imsize=imsize,
          cell=pixsize,
          weighting='natural',
          uvtaper=['200klambda'],
          niter=0)
    dirty_200klam_rms = imstat(imagename='{}.residual'.format(image),
                region='aumic_rms.region', listit=False)['rms'][0]
    if view:
        viewer(infile=image + '.residual')
        print 'dirty taper: {}'.format(dirty_natural_rms)

    # User mask:
    image = filename + "_200klam"
    subprocess.call("rm -rf {}.*".format(image), shell=True)
    tclean(vis=filename + ".concat.ms",
          imagename=image,
          imsize=imsize,
          cell=pixsize,
          weighting='natural',
          uvtaper=['200klambda'],
          niter=100000,
          threshold=dirty_200klam_rms / 2.,
          usemask='user',
          mask=taper_mask,
          pbmask=None)
    if view:
        viewer(infile=image + '.image')
    taper_rms = imstat(imagename='{}.image'.format(image), region='aumic_rms.region',
                listit=False)['rms'][0]
    print 'Taper clean rms is', taper_rms

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
              fitsimage='{}.fits'.format(image))
