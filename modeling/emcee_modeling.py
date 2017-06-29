from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import subprocess as sp
import emcee
cgdisp_start = False

class Observation:
    def __init__(self, name):
        self.name = name
        self.uvf  = fits.open('obs_data/{}.uvf'.format(name))
        
        self.dec = self.uvf[3].data['DECEPO'][0]
        self.ra = self.uvf[3].data['RAEPO'][0]

class Model:
    def obs_sample(self, obs):
        """
        Create model fits file with correct header information and sample using 
        ALMA observation uv coverage to to create model .vis and .uvf files.
        """
        
        print('')
        print('================================================================================')
        print('                                 obs_sample({})                               ').format(obs.name)
        print('================================================================================')
        
        # make observation-specific model name
        filename = self.name + '_' + obs.name

        # First, delete all model files with same name
        sp.call('rm -rf model_data/{}*'.format(filename), shell=True)
        
        # Align with observation and save observation-specific fits file
        self.im[0].header['CRVAL1'] = obs.ra
        self.im[0].header['CRVAL2'] = obs.dec
        self.im.writeto('model_data/{}.fits'.format(filename))

        # Convert model into MIRIAD .im image file
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
        Return chi^2 statistics of model.
        """
        
        print('')
        print('================================================================================')
        print('                                  get_chi({})                                  ').format(obs.name)
        print('================================================================================')

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
        print(np.max(weights))
        chi = np.sum((datrlI - modrlI)**2 * weights +
                     (datimI - modimI)**2 * weights)
        redchi = chi / len(datrlI)

        self.redchis.append(redchi)

    def clean(self, obs, residual=False, show=True):
        """
        Clean and image (if desired) a observation-specific model.
        Either model image or residuals may be chosen.
        """
        
        print('')
        print('================================================================================')
        if residual:
            print('                              clean({}.residual)                               ').format(obs.name)
        else:
            print('                                   clean({})                                   ').format(obs.name)
        print('================================================================================')
        
        # Set observation-specific clean filename; clear filenames
        filename = self.name + '_' + obs.name
        if residual == True:
            filename += '.residual'
        sp.call('rm -rf model_data/{}.{{mp,bm,cl,cm}}'.format(filename), shell=True)
        
        #Dirty clean; save rms for clean cutoff
        sp.call(['invert', 
            'vis=model_data/{}.vis'.format(filename), 
            'map=model_data/{}.mp'.format(filename), 
            'beam=model_data/{}.bm'.format(filename), 
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
        imstat_out=sp.check_output(['imstat', 
            'in=model_data/{}.mp'.format(filename), 
            "region='boxes(256,0,512,200)'"])
        dirty_rms = float(imstat_out[-38:-29])
        print("Dirty rms is {}".format(dirty_rms))
            
        
        # Clean down to half the rms
        sp.call(['clean', 
            'map=model_data/{}.mp'.format(filename), 
            'beam=model_data/{}.bm'.format(filename), 
            'out=model_data/{}.cl'.format(filename), 
            'niters=10000', 'cutoff={}'.format(dirty_rms/2)])
        sp.call(['restor',
            'map=model_data/{}.mp'.format(filename),
            'beam=model_data/{}.bm'.format(filename),
            'model=model_data/{}.cl'.format(filename),
            'out=model_data/{}.cm'.format(filename)])
        
        # Display clean image with 2,4,6 sigma contours, if desired
        if show == True:
            
            # Display an unimportant image to get around the fact that the first 
            # image displayed with cgdisp in a session can't be deleted
            global cgdisp_start
            if cgdisp_start == False:
                sp.call(['cgdisp', 'in=cgdisp_start.im', 'type=p', 'device=/xs'])
                cgdisp_start = True
            
            #Get rms for countours
            imstat_out = sp.check_output(['imstat', 
                'in=model_data/{}.cm'.format(filename), 
                "region='boxes(256,0,512,200)'"])
            clean_rms = float(imstat_out[-38:-29])
            print("Clean rms is {}".format(clean_rms))
            
            # Display
            sp.call(['cgdisp', 
                'in=model_data/{}.cm,model_data/{}.cm'.format(filename, filename), 
                'type=p,c', 'device=/xs', 
                'slev=a,{}'.format(clean_rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])

    def residuals(self, obs, show=True):
        """
        Create model residuals (data - model), and clean//display if desired
        """
        
        #Set observation-specific filename
        filename = self.name + '_' + obs.name
        
        # Subtract model visibilities from data; outfile is residual visibilities
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
        
        # Associate observations with data
        self.observations = observations
        
        # For each observation, append reduced chi^2 to list
        self.redchis = []
        for obs in observations:
            self.obs_sample(obs)
            self.get_chi(obs)
        self.im.close
        
mar = Observation('mar')
aug = Observation('aug')
jun = Observation('jun')
first_model = Model('first_model', [mar,aug,jun])

# Display clean images and residuals for each observation
for obs in first_model.observations:
    first_model.clean(obs)
    raw_input('cool?')
    first_model.residuals(obs)
    raw_input('cool?')
    
# chi = get_chi(aug, model)
