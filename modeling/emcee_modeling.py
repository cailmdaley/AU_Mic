from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess as sp
import emcee

    
class Observation:
    def __init__(self, name):
        self.name = name
        self.im = fits.open('obs_data/{}.fits'.format(name))
        self.uvf  = fits.open('obs_data/{}.uvf'.format(name))

        self.ra = self.im[0].header['CRVAL1']
        self.dec = self.im[0].header['CRVAL2']

class Model:
    
    def obs_sample(self, obs):
        """
        Create model fits file with correct header information, remove stellar emission, convolve with ALMA visiblities to create model .vis and .uvf files.
        """
        
        filename = self.name + '_' + obs.name

        # first, delete all model files with same name
        sp.call('rm -rf model_data/{}*'.format(filename), shell=True)
        
        # align with observation
        self.im[0].header['CRVAL1'] = obs.ra
        self.im[0].header['CRVAL2'] = obs.dec

        self.im.writeto('model_data/{}.fits'.format(filename))

        # Convert model into MIRIAD .im
        sp.call(['fits', 'op=xyin', 
            'in=model_data/{}.fits'.format(filename),
            'out=model_data/{}.im'.format(filename)])

        # Sample the model image using the observation uv coverage
        sp.call(['uvmodel', 'options=replace',
            'vis=obs_data/{}.vis'.format(obs.name), 
            'model=model_data/{}.im'.format(filename),
            'out=model_data/{}.vis'.format(filename)])
            
        #Convert to UVfits
        print(filename)
        sp.call(['fits', 'op=uvout', 
            'in=model_data/{}.vis'.format(filename),
            'out=model_data/{}.uvf'.format(filename)])
            
            
    def get_chi(self, obs):
        """
        Return chi^2 of model.
        """

        datrlimwt = obs.uvf[0].data['data']
        datrlxx = datrlimwt[:, 0, 0, 0, 0, 0, 0]
        datrlyy = datrlimwt[:, 0, 0, 0, 0, 1, 0]
        datimxx = datrlimwt[:, 0, 0, 0, 0, 0, 1]
        datimyy = datrlimwt[:, 0, 0, 0, 0, 1, 1]
        weights = datrlimwt[:, 0, 0, 0, 0, 0, 2]
        datrlI = np.array((datrlxx + datrlyy) / 2.)
        datimI = np.array((datimxx + datimyy) / 2.)

        self.uvf  = fits.open('model_data/{}_{}.uvf'.format(self.name, obs.name))
        modrlimwt = self.uvf[0].data['data']
        modrlI = modrlimwt[::2, 0, 0, 0, 0, 0]
        modimI = modrlimwt[::2, 0, 0, 0, 0, 1]
        self.uvf.close

        # Calculate chi^2
        chi = np.sum((datrlI - modrlI)**2 * weights +
                     (datimI - modimI)**2 * weights)
        redchi = chi / len(datrlI)

        return redchi
            
    def clean(self, obs, residual=False, show=True):
        """
        Clean and image a miriad visibility file; uses imstat to print rms, and then asks the user to input a clean cutoff level.
        """
        filename = self.name + '_' + obs.name
        if residual == True:
            filename += '.residual'

        sp.call('rm -rf model_data/{}.{{mp,bm,cl,cm}}'.format(filename), shell=True)
        sp.call(['invert', 
            'vis=model_data/{}.vis'.format(filename), 
            'map=model_data/{}.mp'.format(filename), 
            'beam=model_data/{}.bm'.format(filename), 
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
            
        sp.call(['imstat', 
            'in=model_data/{}.mp'.format(filename), 
            "region='boxes(256,0,512,200)'"])
        cutoff = input('Please enter CLEAN cutoff:  ')
        # cutoff = 6.657E-12
        
        sp.call(['clean', 
            'map=model_data/{}.mp'.format(filename), 
            'beam=model_data/{}.bm'.format(filename), 
            'out=model_data/{}.cl'.format(filename), 
            'niters=10000', 'cutoff={}'.format(cutoff)])
             
        sp.call(['restor',
            'map=model_data/{}.mp'.format(filename),
            'beam=model_data/{}.bm'.format(filename),
            'model=model_data/{}.cl'.format(filename),
            'out=model_data/{}.cm'.format(filename)])
        
        if show == True:
            
            sp.call(['imstat', 
                'in=model_data/{}.cm'.format(filename), 
                "region='boxes(256,0,512,200)'"])
            rms = input('Please enter CLEAN cutoff:  ')
            
            sp.call(['cgdisp', 
                'in=model_data/{}.cm,model_data/{}.cm'.format(filename, filename), 
                'type=p,c', 'device=/xs', 
                'slev=a,{}'.format(rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])

        
    def residuals(self, obs, show=True):
        filename = self.name + '_' + obs.name
        
        sp.call(['uvmodel', 'options=subtract',
            'model=model_data/{}.im'.format(filename),
            'vis=obs_data/{}.vis'.format(obs.name),
            'out=model_data/{}.residual.vis'.format(filename)])
    
        if show == True:
            self.clean(obs, residual=True)
        
            
    def __init__(self, path, observations, name='model'):
        self.path = path
        self.name = name
        self.im = fits.open(path + '.fits')
        
        self.observations = observations
        
        self.redchis = []
        for obs in observations:
            self.obs_sample(obs)
            self.redchis.append(self.get_chi(obs))
        self.im.close
        

aug = Observation('aug')
jun = Observation('jun')
first_model = Model('first_model', [aug,jun])
# chi = get_chi(aug, model)
