#-------------------------------------------------------------------
# Modeling_Code.py : Modeling Code for ALMA Observations of AU_Mic |
# AUTHOR           : Cail Daley                                    |
# PI               : A. M. Hughes                                  |
# INSTITUTION      : Department of Astronomy,                      |
#                  : Wesleyan University                           |
# DESCRIPTION      : The purpose of this code is...ToDo            |
#===================================================================
from radmc3dPy import *
from astropy.io import fits
import subprocess
import numpy as np
import matplotlib.pylab as plb
#===================================================================#


# Variables and input files:
modelname = 'aumic_model'
filenames = ['18aug2015_aumic_spw0.corrected_weights',
             '18aug2015_aumic_spw1.corrected_weights',
             '18aug2015_aumic_spw2.corrected_weights',
             '18aug2015_aumic_spw3.corrected_weights',
             '24jun2015_aumic1_spw0.corrected_weights',
             '24jun2015_aumic1_spw1.corrected_weights',
             '24jun2015_aumic1_spw2.corrected_weights',
             '24jun2015_aumic1_spw3.timeflag.corrected_weights',
             '26mar2014_aumic_spw0.corrected_weights',
             '26mar2014_aumic_spw1.corrected_weights',
             '26mar2014_aumic_spw2.corrected_weights',
             '26mar2014_aumic_spw3.corrected_weights']
coord = ['20h45m09.854710s -031d20m32.52034s',
         '20h45m09.854710s -031d20m32.52034s',
         '20h45m09.854710s -031d20m32.52034s',
         '20h45m09.854710s -031d20m32.52034s',
         '20h45m09.867700s -031d20m32.89000s',
         '20h45m09.867700s -031d20m32.89000s',
         '20h45m09.867700s -031d20m32.89000s',
         '20h45m09.867700s -031d20m32.89000s',
         '20h45m09.844300s -031d20m31.36000s',
         '20h45m09.844300s -031d20m31.36000s',
         '20h45m09.844300s -031d20m31.36000s',
         '20h45m09.844300s -031d20m31.36000s']
ra = [15. * (20. + 45. / 60. + 9.85471 / 3600.),
      15. * (20. + 45. / 60. + 9.85471 / 3600.),
      15. * (20. + 45. / 60. + 9.85471 / 3600.),
      15. * (20. + 45. / 60. + 9.85471 / 3600.),
      15. * (20. + 45. / 60. + 9.867700 / 3600.),
      15. * (20. + 45. / 60. + 9.867700 / 3600.),
      15. * (20. + 45. / 60. + 9.867700 / 3600.),
      15. * (20. + 45. / 60. + 9.867700 / 3600.),
      15. * (20. + 45. / 60. + 9.844300 / 3600.),
      15. * (20. + 45. / 60. + 9.844300 / 3600.),
      15. * (20. + 45. / 60. + 9.844300 / 3600.),
      15. * (20. + 45. / 60. + 9.844300 / 3600.)]
dec = [-31. - 20. / 60. - 32.52034 / 3600.,
       -31. - 20. / 60. - 32.52034 / 3600.,
       -31. - 20. / 60. - 32.52034 / 3600.,
       -31. - 20. / 60. - 32.52034 / 3600.,
       -31. - 20. / 60. - 32.89 / 3600.,
       -31. - 20. / 60. - 32.89 / 3600.,
       -31. - 20. / 60. - 32.89 / 3600.,
       -31. - 20. / 60. - 32.89 / 3600.,
       -31. - 20. / 60. - 32.36 / 3600.,
       -31. - 20. / 60. - 32.36 / 3600.,
       -31. - 20. / 60. - 32.36 / 3600.,
       -31. - 20. / 60. - 32.36 / 3600.]
starflux = 9.417385e-05
pixsize = '0.03arcsec'


def create_model():
    """
    Returns a model image of a circumstellar disk.
    """

    analyze.writeDefaultParfile('ppdisk')
    setup.problemSetupDust('ppdisk', crd_sys='sph', nx=[20], ny=[20], nz=[20], xbound='[0.1*au,100*au]', ybound='[0.0,pi]', zbound='[0.0,2*pi]', nw=[10], wbound=['0.1', '1500'], bgdens=1e-30, dusttogas=0.01, gasspec_vturb=1000, gsmax='10', gsmin='1',
                           hpr_prim_rout=0.0, hrdisk=0.1, hrpivot='100*au', itempdecoup=0, lines_mode=0, mdisk=1e28, mstar='0.31*ms', nphot=1e6, plh=1 / .7, plsig1=2.8, rin='8.8*au', rdisk='40.3*au', rstar='0.86*rs', sig0=0.0, srim_plsig=0.0, srim_rout=0.0, tstar=3600.)

    par = analyze.readParams()
    par.printPar()

    subprocess.call('radmc3d mctherm', shell=True)
    image.makeImage(npix=300., wav=1351.0, incl=89.5,
                    posang=218.41, sizeau=150.)
    im = image.readImage()

    return im


def model_convolve(im, modelname, filename, coord, ra, dec):
    """
    Create model fits file with correct header information, remove stellar emission, convolve with ALMA visiblities to create model .vis and .uvf files.
    """

    subprocess.call('rm -r model_files', shell=True)
    subprocess.call('mkdir model_files', shell=True)

    im.writeFits('model_files/{}.fits'.format(modelname),
                 dpc=8.9, coord='{}'.format(coord))
    model = fits.open('model_files/{}.fits'.format(modelname))
    model[0].data[0, 150, 150] -= starflux
    model[0].header['CRVAL1'] = ra
    model[0].header['CRVAL2'] = dec

    model.writeto('model_files/{}.fits'.format(modelname), overwrite=True)
    model.close
    subprocess.call(['fits', 'in=model_files/{}.fits'.format(modelname),
                     'op=xyin', 'out=model_files/{}.im'.format(modelname)])
    subprocess.call(['uvmodel', 'model=model_files/{}.im'.format(modelname), 'vis=data_files/{}.vis'.format(
        filename), 'options=replace', 'out=model_files/{}.vis'.format(modelname)])
    subprocess.call(['fits', 'in=model_files/{}.vis'.format(modelname),
                     'op=uvout', 'out=model_files/{}.uvf'.format(modelname)])


def get_chi(filename, modelname):
    """
    Return chi^2 of model.
    """

    data = fits.open('data_files/{}.uvf'.format(filename))
    datrlimwt = data[0].data['data']

    # splitting visibilities by time removes one index from datrlimwt array:
    # so, must check size of datrlimwt
    if len(datrlimwt.shape) == 7:
        datrlxx = datrlimwt[:, 0, 0, 0, 0, 0, 0]
        datrlyy = datrlimwt[:, 0, 0, 0, 0, 1, 0]
        datimxx = datrlimwt[:, 0, 0, 0, 0, 0, 1]
        datimyy = datrlimwt[:, 0, 0, 0, 0, 1, 1]
        weights = datrlimwt[:, 0, 0, 0, 0, 0, 2]
    elif len(datrlimwt.shape) == 6:
        datrlxx = datrlimwt[:, 0, 0, 0, 0, 0]
        datrlyy = datrlimwt[:, 0, 0, 0, 1, 0]
        datimxx = datrlimwt[:, 0, 0, 0, 0, 1]
        datimyy = datrlimwt[:, 0, 0, 0, 1, 1]
        weights = datrlimwt[:, 0, 0, 0, 0, 2]
    datrlI = np.array((datrlxx + datrlyy) / 2.)
    datimI = np.array((datimxx + datimyy) / 2.)

    model = fits.open('model_files/{}.uvf'.format(modelname))
    modrlimwt = model[0].data['data']
    modrlI = modrlimwt[::2, 0, 0, 0, 0, 0]
    modimI = modrlimwt[::2, 0, 0, 0, 0, 1]

    data.close
    model.close


# Calculate chi^2
    chi = np.sum((datrlI - modrlI)**2 * weights +
                 (datimI - modimI)**2 * weights)
    redchi = chi / len(datrlI)

    return chi, redchi


def image_vis(vis, show=True):
    """
    Clean and image a miriad visibility file; uses imstat to print rms, and then asks the user to input a clean cutoff level.
    """

    subprocess.call(['rm -r ' + vis + '.{mp,bm,cl,cm}'], shell=True)
    subprocess.call('invert vis={}.vis map={}.mp beam={}.bm cell={} imsize=512 options=systemp,mfs robust=2'.format(
        vis, vis, vis, pixsize), shell=True)
    subprocess.call(['imstat', 'in={}.mp'.format(vis),
                     "region='boxes(256,0,512,200)'"])
    cutoff = input('Please enter CLEAN cutoff:  ')
    subprocess.call('clean map={}.mp beam={}.bm out={}.cl niters=10000 cutoff={}'.format(
        vis, vis, vis, cutoff), shell=True)  # cutoff=6.290E-05
    subprocess.call('restor map={}.mp beam={}.bm model={}.cl out={}.cm'.format(
        vis, vis, vis, vis), shell=True)
    if show == True:
        subprocess.call(['cgdisp', 'in={}.cm'.format(
            vis), 'device=/xs', 'labtyp=arcsec', 'beamtyp=b,l,3'])


im = create_model()

# Calculate chi^2 of each model
chis = []
redchis = []
for i in range(len(filenames)):
    model_convolve(im, modelname + str(i),
                   filenames[i], coord[11], ra[11], dec[11])
    chi, redchi = get_chi(filenames[i], modelname + str(i))

    chis.append(chi)
    redchis.append(redchi)
