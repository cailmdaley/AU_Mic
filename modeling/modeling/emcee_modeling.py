from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import emcee


def model_convolve(modelfile, obsfile, coord, ra, dec):
    """
    Create model fits file with correct header information, remove stellar emission, convolve with ALMA visiblities to create model .vis and .uvf files.
    """

    # subprocess.call('rm -r model_files', shell=True)
    # subprocess.call('mkdir model_files', shell=True)

    #Read in model and observations
    # model = fits.open('{}.fits'.format(modelname))[0]
    model = fits.open('first_model.fits')[0]
    obs = fits.open('{}.fits'.format(obsfile))[0]
    
    #Align pointing centers of model and obs
    model[0].header['CRVAL1'] = obs[0].header['CRVAL1']
    model[0].header['CRVAL2'] = obs[0].header['CRVAL2']
    
    #Write aligned model
    model.writeto('model_files/{}.fits'.format(modelname), overwrite=True)
    model.close
        
    subprocess.call(['fits', 'in=model_files/{}.fits'.format(modelname),
                     'op=xyin', 'out=model_files/{}.im'.format(modelname)])
    subprocess.call(['uvmodel', 'model=model_files/{}.im'.format(modelname), 
        'vis=../data_files/{}.vis'.format(datafile), 'options=replace', 
        'out=model_files/{}.vis'.format(modelname)])
    subprocess.call(['fits', 'in=model_files/{}.vis'.format(modelname),
        'op=uvout', 'out=model_files/{}.uvf'.format(modelname)])

                     
# model_data = fits.getdata('first_model.fits')
# model_head = fits.getheader('first_model.fits')

blah = fits.open('first_model.fits')[0]
blah.data.shape
blah[0].data[0]

filename = 'first_model.fits'
filename[:-5]
obs_data = fits.open('aumic_all_centered_natural.fits')[0]
obs_data.header

def chi(model_data):
    return (obs_data - model_data)**2 / rms

model_data[0].shape
plt.imshow(model_data[0], origin='lower')
plt.show()
