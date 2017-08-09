import numpy as np
import subprocess as sp
from astropy.io import fits
from collections import OrderedDict

from astrocail import fitting, plotting
from disk_model import debris_disk, raytrace
from aumic_observations import band6_observations, band6_rms_values, band6_fits_images


# default parameter dict:
param_dict = OrderedDict([
    ('temp_index',        -0.5),
    ('m_disk',            -8.),
    ('sb_law',            2.3),
    ('r_in',              8.8),
    ('d_r',                31.5),
    ('r_crit',            150.0),
    ('inc',               89.5),
    ('m_star',            0.31),
    ('co_frac',           0.0001),
    ('v_turb',            0.081),
    ('Zq',                70.0),
    ('column_densities', [0.79, 1000]),
    ('abundance_bounds', [50, 500]),
    ('hand',              -1),
    ('rgrid_size',        500),
    ('zgrid_size',        500),
    ('l_star',            0.09),
    ('scale_factor',      0.1),
    ('pa',                128.41),
    ('starflux',      2.50e-4)])

def make_fits(model, disk_params):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(structure_params, obs=[300, 131, 300, 20])
    raytrace.total_model(model_disk,
        distance=9.91, # pc
        imres=0.03, # arcsec/pix
        xnpix=512, #image size in pixels
        freq0=model.observations[0][0].uvf[0].header['CRVAL4']*1e-9, # obs frequeency
        PA=PA,
        offs=[0.0,0.0], # offset from image center
        nchans=1, # continum
        isgas=False, # continuum!
        includeDust=True, #continuuum!!
        extra=0.0, # ?
        modfile = model.root + model.name)

    
def fix_fits(model, obs, starflux):
    # open fits and add starflux,
    model_fits = fits.open(model.path + '.fits')
    crpix = int(model_fits[0].header['CRPIX1'])
    model_im = model_fits[0].data[0]
    try:
        model_im[crpix, crpix] = model.crpix_diskflux + starflux
    except AttributeError:
        model.crpix_diskflux = model_im[crpix, crpix]
        model_im[crpix, crpix] += starflux
    
    
    model_fits[0].header['CRVAL1'] = obs.ra
    model_fits[0].header['CRVAL2'] = obs.dec
    
    model_fits.writeto(model.path + '.fits', overwrite=True)
        
# define likelehood functions
def lnprob(theta, run_name, to_vary):
    """
    For each parameter to be varied, return -infinity if the proposed value lies
    outside the bounds for the parameter; if all parameters are ok, return 0.
    """
    
    # reassign parameters to vary to those input by emcee
    for i, free_param in enumerate(to_vary):
        lower_bound, upper_bound = free_param[-1]
        
        if lower_bound < theta[i] < upper_bound:
            param_dict[free_param[0]] = theta[i]
        else: return -np.inf

    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']
    
    disk_params = param_dict.values()[:-1]
    starflux = param_dict.values()[-1]
    
    # intialize model and make fits image 
    model = fitting.Model(observations=band6_observations,
        root=run_name + '/model_files/', 
        name='model' + str(np.random.randint(1e10)))
    make_fits(model, disk_params)
    for obs in band6_observations.flatten():
        fix_fits(model, obs, starflux)
        model.obs_sample(obs)
        model.get_chi(obs)
    
    model.delete()
    return -0.5 * sum(model.chis)


def make_best_fits(run):
    best_fit = run.chain[run.chain['lnprob'] == run.chain['lnprob'].max()]
    print('Best-fit model parameters:')
    print(best_fit.to_string())
    print('')
        
    for param in best_fit.columns[:-1]:
        param_dict[param] = best_fit[param].values
        
        
    param_dict['m_disk'] = 10**param_dict['m_disk']
    param_dict['d_r'] += param_dict['r_in']

    disk_params = param_dict.values()[:-1]
    starflux = param_dict.values()[-1]
    
    # intialize model and make fits image 
    print('Making model...')
    model = fitting.Model(observations=band6_observations,
        root=run.name + '/model_files/', 
        name=run.name + '_bestfit')
    # model.delete()
    # make_fits(model, disk_params)
    
    print('Sampling...')
    visibilities = []
    residuals = []
    for pointing, rms in zip(band6_observations, band6_rms_values):
        ids = []
        for obs in pointing:
            pass
            # fix_fits(model, obs, starflux)
            # 
            # ids.append('_' + obs.name[12:20])
            # model.obs_sample(obs, ids[-1])
            # model.make_residuals(obs, ids[-1])
            
        # make model visibilities
        # vis_files = ','.join([model.path+ident+'.vis' for ident in ids])
        visibilities.append('{}_{}'.format(model.path, obs.name[12:15]))
        # sp.call(['uvcat', 'vis={}'.format(vis_files), 'out={}.vis'.format(visibilities[-1])])
        # model.clean(visibilities[-1], rms, show=False)
        
        # make residuals
        # res_files = ','.join([model.path+ident+'.residuals.vis' for ident in ids])
        residuals.append('{}_{}.residuals'.format(model.path, obs.name[12:15]))
        # sp.call(['uvcat', 'vis={}'.format(res_files), 'out={}.vis'.format(residuals[-1])])
        # model.clean(residuals[-1], rms, show=False)
        
    fig = plotting.Figure(
        [[band6_fits_images[i], visibilities[i]+'.fits', residuals[i]+'.fits'] for i in range(len(visibilities))],
        [3*[rms] for rms in band6_rms_values])
        
        
        
    
        
    
if __name__ == '__main__':
    fitting.run_emcee(run_name='run7', nsteps=10000, nwalkers=18, 
        lnprob = lnprob, to_vary = [
        ('m_disk',             -7.55,         0.05,       (-np.inf, np.inf)),
        ('sb_law',             2.3,           2,          (-5.,     10.)), 
        ('scale_factor',       0.05,          0.03,       (0,       np.inf)),
        ('r_in',               8.8,           5,          (0,       np.inf)),
        ('d_r',                31.5,          10,         (0,       np.inf)),
        ('inc',                90,            1,          (0,       np.inf)),
        ('pa',                 128.48,        0.1,        (0,       360)),
        ('starflux',           2.50e-4,       1e-4,      (0,       np.inf))])
