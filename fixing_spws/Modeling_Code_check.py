# Modeling Code
from radmc3dPy import *
from astropy.io import fits
import subprocess
import numpy as np
import matplotlib.pylab as plb
#==============================================================================#


# Variables and input files:
modelname = 'aumic_model'
filenames = ['18aug2015_aumic_spw0.corrected_weights',
             '18aug2015_aumic_spw1.corrected_weights',
             '18aug2015_aumic_spw2.corrected_weights',
             '18aug2015_aumic_spw3.corrected_weights',
             '24jun2015_aumic1_spw0.corrected_weights.timesplit',
             '24jun2015_aumic1_spw1.corrected_weights.timesplit',
             '24jun2015_aumic1_spw2.corrected_weights.timesplit',
             '24jun2015_aumic1_spw3.corrected_weights.timesplit',
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

# Create Model


def create_model():
    """
    Returns a model image of the disk
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

# Add star position & subtract star flux, create convolved model visibilities


def model_convolve(im, modelname, filename, coord, ra, dec):
    """
    Create model fits file with correct header information, remove stellar emission, convolve with ALMA visiblities to create model .vis and .uvf files.
    """
    print(modelname)

    subprocess.call(['rm -rf {}*'.format(modelname)], shell=True)

    im.writeFits('{}.fits'.format(modelname),
                 dpc=8.9, coord='{}'.format(coord))
    model = fits.open('{}.fits'.format(modelname))
    model[0].data[0, 150, 150] -= starflux
    model[0].header['CRVAL1'] = ra
    model[0].header['CRVAL2'] = dec

    # Clear pre-existing models before writing
    subprocess.call(['rm -rf {}.fits'.format(modelname)], shell=True)
    model.writeto('{}.fits'.format(modelname))
    model.close
    subprocess.call(['fits', 'in={}.fits'.format(modelname),
                     'op=xyin', 'out={}.im'.format(modelname)])
    subprocess.call(['uvmodel', 'model={}.im'.format(modelname), 'vis={}.vis'.format(
        filename), 'options=replace', 'out={}.vis'.format(modelname)])
    subprocess.call(['fits', 'in={}.vis'.format(modelname),
                     'op=uvout', 'out={}.uvf'.format(modelname)])


def get_chi(filename, modelname):
    """
    Return chi^2 of model using seperate weightfile.
    """

    data = fits.open('{}.uvf'.format(filename))
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

    model = fits.open('{}.uvf'.format(modelname))
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

# Make this part do the file splitting
#subprocess.call("uvlist vis=24jun2015_aumic1_spw3.vis options=baseline")
#os.system("prthd in=24jun2015_aumic1_spw3.vis")


def vis_cut(vis, cut, suff, filenames):
    subprocess.call('rm -r {}.vis'.format(vis + suff), shell=True)
    subprocess.call(['uvaver vis={}.vis select={} out={}.vis'.format(
        vis, cut, vis + suff)], shell=True)
    filenames.append(vis + suff)


def create_uvf(filename):
    subprocess.call('rm -r {}.uvf'.format(filename), shell=True)
    subprocess.call('fits in={}.vis op=uvout out={}.uvf'.format(
        filename, filename), shell=True)


def image_vis(vis, pixsize, show=True):
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


visibilities = ['24jun2015_aumic1_spw0', '24jun2015_aumic1_spw1', '24jun2015_aumic1_spw2', '24jun2015_aumic1_spw3']
newfiles = []
for vis in visibilities:
    vis_cut(vis, "'time(15JUN24:03:45:36.0,15JUN24:04:20:00.0)'",
            ".timesplit", newfiles)
    create_uvf(vis+'.timesplit')

im = create_model()

chis = []
redchis = []
for i in range(len(filenames)):
    model_convolve(im, modelname + str(i),
                   filenames[i], coord[11], ra[11], dec[11])
    chi, redchi = get_chi(filenames[i], modelname + str(i))

    chis.append(chi)
    redchis.append(redchi)
