import emcee
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
from disk_model import debris_disk
from astropy.io import fits
cgdisp_start = False

class Observation:
    def __init__(self, name, rms):
        self.name = name
        self.uvf  = fits.open('obs_data/{}.uvf'.format(name))
        
        self.rms = rms
        
        self.dec = self.uvf[3].data['DECEPO'][0]
        self.ra = self.uvf[3].data['RAEPO'][0]
        
    def clean(self, show=True):
        """
        Clean and image (if desired) a observation-specific model.
        Either model image or residuals may be chosen.
        """
        
        print('')
        print('================================================================================')
        print('                            clean({})                                   ').format(self.name)
        print('================================================================================')
        
        # Set observation-specific clean filename; clear filenames
        sp.call('rm -rf {}.{{mp,bm,cl,cm}}'.format(self.name), shell=True)
        
        #Dirty clean; save rms for clean cutoff
        sp.call(['invert', 
            'vis=obs_data/{}.vis'.format(self.name), 
            'map=obs_data/{}.mp'.format(self.name), 
            'beam=obs_data/{}.bm'.format(self.name), 
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
        imstat_out=sp.check_output(['imstat', 
            'in=obs_data/{}.mp'.format(self.name), 
            "region='boxes(256,0,512,200)'"])
        dirty_rms = float(imstat_out[-38:-29])
        print("Dirty rms is {}".format(dirty_rms))
            
        
        # Clean down to half the rms
        sp.call(['clean', 
            'map=obs_data/{}.mp'.format(self.name), 
            'beam=obs_data/{}.bm'.format(self.name), 
            'out=obs_data/{}.cl'.format(self.name), 
            'niters=10000', 'cutoff={}'.format(dirty_rms/2)])
        sp.call(['restor',
            'map=obs_data/{}.mp'.format(self.name),
            'beam=obs_data/{}.bm'.format(self.name),
            'model=obs_data/{}.cl'.format(self.name),
            'out=obs_data/{}.cm'.format(self.name)])
        
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
                'in=obs_data/{}.cm'.format(self.name), 
                "region='boxes(256,0,512,200)'"])
            clean_rms = float(imstat_out[-38:-29])
            print("Clean rms is {}".format(clean_rms))
            
            # Display
            sp.call(['cgdisp', 
                'in=obs_data/{}.cm,obs_data/{}.cm'.format(self.name, self.name), 
                'type=p,c', 'device=/xs', 
                'slev=a,{}'.format(clean_rms), 'levs1=-6,-4,-2,2,4,6',
                'region=arcsec,box(-5,-5,5,5)',
                'labtyp=arcsec', 'beamtyp=b,l,3',])

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
        sp.call(['fits', 'op=uvout', 
            'in=model_data/{}.vis'.format(filename),
            'out=model_data/{}.uvf'.format(filename)])
            
    def get_chi(self, obs):
        """
        Return chi^2 statistics of model.
        """
        
        datrlimwt = obs.uvf[0].data['data']
        datrl_xx = datrlimwt[:, 0, 0, 0, 0, 0, 0]
        datrl_yy = datrlimwt[:, 0, 0, 0, 0, 1, 0]
        datim_xx = datrlimwt[:, 0, 0, 0, 0, 0, 1]
        datim_yy = datrlimwt[:, 0, 0, 0, 0, 1, 1]
        weights =  datrlimwt[:, 0, 0, 0, 0, 0, 2]
        datrl_stokes = np.array((datrl_xx + datrl_yy) / 2.)
        datim_stokes = np.array((datim_xx + datim_yy) / 2.)

        self.uvf  = fits.open('model_data/{}_{}.uvf'.format(self.name, obs.name))
        modrlimwt = self.uvf[0].data['data']
        modrl_stokes = modrlimwt[::2, 0, 0, 0, 0, 0]
        modim_stokes = modrlimwt[::2, 0, 0, 0, 0, 1]
        self.uvf.close

        # Calculate chi^2
        chi = np.sum((datrl_stokes - modrl_stokes)**2 * weights +
                     (datim_stokes - modim_stokes)**2 * weights)
        redchi = chi / len(datrl_stokes)

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
        
        # Clean down to half the observation rms
        sp.call(['invert', 
            'vis=model_data/{}.vis'.format(filename), 
            'map=model_data/{}.mp'.format(filename), 
            'beam=model_data/{}.bm'.format(filename), 
            'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
        sp.call(['clean', 
            'map=model_data/{}.mp'.format(filename), 
            'beam=model_data/{}.bm'.format(filename), 
            'out=model_data/{}.cl'.format(filename), 
            'niters=10000', 'cutoff={}'.format(obs.rms/2)])
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
        
mar0 = Observation('aumic_mar_spw0_FINAL', rms=6.5e-05)
mar1 = Observation('aumic_mar_spw1_FINAL', rms=6.124e-05)
mar2 = Observation('aumic_mar_spw2_FINAL', rms=6.068e-05)
mar3 = Observation('aumic_mar_spw3_FINAL', rms=6.468e-05)
aug0 = Observation('aumic_aug_spw0_FINAL', rms=5.879e-05)
aug1 = Observation('aumic_aug_spw1_FINAL', rms=5.336e-05)
aug2 = Observation('aumic_aug_spw2_FINAL', rms=6.092e-05)
aug3 = Observation('aumic_aug_spw3_FINAL', rms=5.558e-05)
jun0 = Observation('aumic_jun_spw0_FINAL', rms=5.369e-05)
jun1 = Observation('aumic_jun_spw1_FINAL', rms=4.658e-05)
jun2 = Observation('aumic_jun_spw2_FINAL', rms=5.083e-05)
jun3 = Observation('aumic_jun_spw3_FINAL', rms=5.559e-05)
observations=[mar0, mar1, mar2, mar3,
              aug0, aug1, aug2, aug3,
              jun0, jun1, jun2, jun3]

first_model = Model('first_model', observations)

# Display clean images and residuals for each observation
# for obs in first_model.observations:
#     first_model.clean(obs)
#     raw_input('cool?')
#     first_model.residuals(obs)
#     raw_input('cool?')
