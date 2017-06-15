from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import emcee

class Observation:
    def __init__(self, name):
        self.name = name
        self.fits = fits.open(name + '.fits')
        self.uvf  = fits.open(name + '.uvf')
        # self.visname  = fits.open(name + '.uvf')

        self.ra = self.fits[0].header['CRVAL1']
        self.dec = self.fits[0].header['CRVAL2']

class Model:
    def __init__(self, name):
        self.name = name
        self.fits = fits.open(name + '.fits')

def chi(model_data):
    return (obs_data - model_data)**2 / rms


def model_convolve(model, obs):

    """
    Create model fits file with correct header information, remove stellar emission, convolve with ALMA visiblities to create model .vis and .uvf files.
    """

    # subprocess.call('mkdir model_files', shell=True)

    #Read in model and observations
    # model = fits.open('{}.fits'.format(modelname))[0]

    #If model lacks a pointing center, align with observation
    if model.fits[0].header['CRVAL1'] == 0:
        model.fits[0].header['CRVAL1'] = obs.ra
        model.fits[0].header['CRVAL2'] = obs.dec

        model.fits.writeto('{}.fits'.format(model.name), clobber=True)
        model.fits.close

    subprocess.call('rm -r {}.im'.format(model.name), shell=True)
    subprocess.call(['fits', 'in={}.fits'.format(model.name),
                     'op=xyin', 'out={}.im'.format(model.name)])

    subprocess.call('rm -r {}.vis'.format(model.name), shell=True)
    subprocess.call(['uvmodel', 'model={}.im'.format(model.name),
        'vis={}.vis'.format(obs.name), 'options=replace',
        'out={}.vis'.format(model.name)])
    subprocess.call(['fits', 'in={}.vis'.format(model.name),
        'op=uvout', 'out={}.uvf'.format(model.name)])

def get_chi(obs, model):
    """
    Return chi^2 of model.
    """

    datrlimwt = obs.uvf[0].data['data']

    # splitting visibilities by time removes one index from datrlimwt array:
    # so, must check size of datrlimwt
    datrlxx = datrlimwt[:, 0, 0, 0, 0, 0, 0]
    datrlyy = datrlimwt[:, 0, 0, 0, 0, 1, 0]
    datimxx = datrlimwt[:, 0, 0, 0, 0, 0, 1]
    datimyy = datrlimwt[:, 0, 0, 0, 0, 1, 1]
    weights = datrlimwt[:, 0, 0, 0, 0, 0, 2]
    # elif len(datrlimwt.shape) == 6:
    #     datrlxx = datrlimwt[:, 0, 0, 0, 0, 0]
    #     datrlyy = datrlimwt[:, 0, 0, 0, 1, 0]
    #     datimxx = datrlimwt[:, 0, 0, 0, 0, 1]
    #     datimyy = datrlimwt[:, 0, 0, 0, 1, 1]
    #     weights = datrlimwt[:, 0, 0, 0, 0, 2]
    datrlI = np.array((datrlxx + datrlyy) / 2.)
    datimI = np.array((datimxx + datimyy) / 2.)

    model.uvf  = fits.open(model.name + '.uvf')
    modrlimwt = model.uvf[0].data['data']
    modrlI = modrlimwt[::2, 0, 0, 0, 0, 0]
    modimI = modrlimwt[::2, 0, 0, 0, 0, 1]

    # Calculate chi^2
    chi = np.sum((datrlI - modrlI)**2 * weights +
                 (datimI - modimI)**2 * weights)
    redchi = chi / len(datrlI)

    print chi, redchi
    # obs.close
    # model.close

def image_vis(model, show=True):
    """
    Clean and image a miriad visibility file; uses imstat to print rms, and then asks the user to input a clean cutoff level.
    """

    pixsize = '0.03arcsec'
    subprocess.call(['rm -r ' + model.name + '.{mp,bm,cl,cm}'], shell=True)
    subprocess.call('invert vis={}.vis map={}.mp beam={}.bm cell={} imsize=512 options=systemp,mfs robust=2'.format(
        model.name, model.name, model.name, pixsize), shell=True)
    subprocess.call(['imstat', 'in={}.mp'.format(model.name),
                     "region='boxes(256,0,512,200)'"])
    cutoff = input('Please enter CLEAN cutoff:  ')
    subprocess.call('clean map={}.mp beam={}.bm out={}.cl niters=10000 cutoff={}'.format(
        model.name, model.name, model.name, cutoff), shell=True)  # cutoff=6.290E-05
    subprocess.call('restor map={}.mp beam={}.bm model={}.cl out={}.cm'.format(
        model.name, model.name, model.name, model.name), shell=True)
    if show == True:
        subprocess.call(['cgdisp', 'in={}.cm'.format(
            model.name), 'device=/xs', 'labtyp=arcsec', 'beamtyp=b,l,3'])


aug = Observation('obs_data/aumic_aug')
model = Model('model_data/first_model')
model_convolve(model, aug)
chi = get_chi(aug, model)
