from astrocail.fitting import *
from astrocail.mcmc import *
from aumic_observations import band6_observations
from disk_model import debris_disk, raytrace
from collections import OrderedDict


# default parameter dict:
params = OrderedDict([
    ('temp_index',        -0.5),
    ('m_disk',            3.67e-08),
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

def make_fits(disk_params, model):
    structure_params = disk_params[:-1]
    PA = disk_params[-1]

    model_disk = debris_disk.Disk(structure_params, obs=[300, 131, 300, 20])
    raytrace.total_model(model_disk,
        distance=9.91, # pc
        imres=0.03, # arcsec/pix
        xnpix=512, #image size in pixels
        freq0=model.observations[0].uvf[0].header['CRVAL4']*1e-9, # obs frequeency
        PA=PA,
        offs=[0.0,0.0], # offset from image center
        nchans=1, # continum
        isgas=False, # continuum!
        includeDust=True, #continuuum!!
        extra=0.0, # ?
        modfile = model.root + model.name)

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
            params[free_param[0]] = theta[i]
        else: return -np.inf

    params['m_disk'] = 10**params['m_disk']
    params['d_r'] += params['r_in']
    
    disk_params = params.values()[:-1]
    starflux = params.values()[-1]
    
    # intialize model and make fits image 
    model = Model(observations=band6_observations,
        root=run_name + '/model_files/', 
        name='model' + str(np.random.randint(1e10)))
    make_fits(disk_params, model)
    
    
    # open fits and add starflux,
    model_fits = fits.open(model.root + model.name + '.fits')
    crpix = int(model_fits[0].header['CRPIX1'])
    model_im = model_fits[0].data[0]
    model_im[crpix, crpix] += starflux
    
    # align with each observation and get chi^2
    for obs in model.observations:
        model_fits[0].header['CRVAL1'] = obs.ra
        model_fits[0].header['CRVAL2'] = obs.dec

        #save and get chi^2
        model_fits.writeto(model.path + '.fits', overwrite=True)
        model.obs_sample(obs)
        model.get_chi(obs)
    
    model.delete()
    return -0.5 * sum(model.chis)

#==============================================================================#
# Sandbox
#==============================================================================#

run_emcee(run_name='run7', nsteps=10000, nwalkers=18, lnprob=lnprob, to_vary = [
    ('m_disk',             -7.55,         0.05,       (-np.inf, np.inf)),
    ('sb_law',             2.3,           2,          (-5.,     10.)), 
    ('scale_factor',       0.05,          0.03,       (0,       np.inf)),
    ('r_in',               8.8,           5,          (0,       np.inf)),
    ('d_r',                31.5,          10,         (0,       np.inf)),
    ('inc',                90,            1,          (0,       np.inf)),
    ('pa',                 128.48,        0.1,        (0,       360)),
    ('starflux',           2.50e-4,       1e-4,      (0,       np.inf))])

# Display clean images and residuals for each observation
# standard_model = Model(params.values(), band6_observations)
# for obs in standard_model.observations:
#     standard_model.clean(obs)
#     raw_input('cool?')
#     standard_model.residuals(obs)
#     raw_input('cool?')
