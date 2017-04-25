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
concat_files = [dates[1], dates[2]]
print 'Files to be cleaned are:', concat_files
raw_input('ok?: ')
# concat_files.append('../Visibility_processing/24jun2015_aumic1_spw3.timeflag.reweighted.ms')

filename = 'aumic_marjun'
imsize = 512
pixsize ='0.03arcsec',

natural = True
natural_mask = 'aumic_marjun_long_bsl_mask.region'

taper = True
taper_mask = 'aumic_marjun_short_bsl_mask.region'

# Concat before clean to fix proper motion offset;
#pointing center of first chronological date is used
if len(concat_files) > 1:
    print('Concatenating...')
    subprocess.call("rm -rf {}".format(filename + ".concat.ms"), shell=True)
    concat(vis=concat_files, concatvis=(filename + ".concat.ms"), dirtol='2arcsec' )

#==========================================================
#rms values for various cleans
aug_natural_rms = 3.02602056763e-05
jun_timing_rms = 2.09834179259e-05
augjun_natural_rms = 1.86735014722e-05
mar_200klam_rms = 2.94901055895e-05
composite_natural_rms = 1.51412368723e-05
composite_200klam_rms = 1.91825256479e-05
marjun_natural_rms = 1.66271911439e-05
marjun_200klam_rms = 2.13612966036e-05
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
    rms = imstat(imagename='{}.residual'.format(image),
                region='aumic_rms.region', listit=False)['rms'][0]

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
          threshold=rms / 2.,
          usemask='user',
          mask=mask,
          pbmask=None)
    viewer(infile=image + '.image')
    taper_rms = imstat(imagename='{}.image'.format(image), region='aumic_rms.region',
                listit=False)['rms'][0]
    print 'Taper clean rms is', taper_rms

    # Export to .fits
    exportfits(imagename='{}.image'.format(image),
              fitsimage='{}.fits'.format(image))
