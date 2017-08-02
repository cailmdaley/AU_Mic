from astrocail.fitting import *
from disk_model import debris_disk, raytrace
from collections import OrderedDict
from emcee.utils import MPIPool
import pandas as pd
import emcee

#==============================================================================#
# Create observations, default parameter dict
#==============================================================================#
mar0 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw0_FINAL', rms=6.5e-05)
mar1 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw1_FINAL', rms=6.124e-05)
mar2 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw2_FINAL', rms=6.068e-05)
mar3 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_mar_spw3_FINAL', rms=6.468e-05)
aug0 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw0_FINAL', rms=5.879e-05)
aug1 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw1_FINAL', rms=5.336e-05)
aug2 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw2_FINAL', rms=6.092e-05)
aug3 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_aug_spw3_FINAL', rms=5.558e-05)
jun0 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw0_FINAL', rms=5.369e-05)
jun1 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw1_FINAL', rms=4.658e-05)
jun2 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw2_FINAL', rms=5.083e-05)
jun3 = Observation(root='unsubtracted_obs_files/', name='aumic_band6_jun_spw3_FINAL', rms=5.559e-05)
band6_observations=[[mar0, mar1, mar2, mar3],
                    [aug0, aug1, aug2, aug3],
                    [jun0, jun1, jun2, jun3]]

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
    ('mar_starflux',      3.67e-4),
    ('aug_starflux',      1.23e-4),
    ('jun_starflux',      2.62e-4)])

def make_fits(disk_params, model):
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

    # create model and get chi^2
    disk_params = params.values()[:-3]
    star_fluxes = params.values()[-3:]

    
    model = Model(observations=band6_observations,
        root=run_name + '/model_files/', 
        name='model' + str(np.random.randint(1e10)))
    make_fits(disk_params, model)

    for date, star_flux in zip(model.observations, star_fluxes):
        for obs in date:
            model_fits = fits.open(model.root + model.name + '.fits')

            # add starflux to central pixel
            crpix = int(model_fits[0].header['CRPIX1'])
            model_im = model_fits[0].data[0]
            model_im[crpix, crpix] += star_flux

            # align with observation
            model_fits[0].header['CRVAL1'] = obs.ra
            model_fits[0].header['CRVAL2'] = obs.dec

            #save
            model_fits.writeto(model.root + model.name + obs.name + '.fits')

            model.obs_sample(obs)
            model.get_chi(obs)

    return -0.5 * sum(model.chis)

def run_emcee(nsteps, nwalkers, run_name, to_vary):
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    # initiate sampler chain
    ndim = len(to_vary)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(run_name, to_vary), pool=pool)

    try:
        samples = pd.read_csv(run_name + '/samples.csv')
        print('Resuming {} at step {}'.format(run_name, samples.index[-1]//nwalkers))
        pos = np.array(samples.iloc[-nwalkers:, :-1])
        with open(run_name + '.csv', 'a') as f:
            f.write('\n')
    except IOError:
        sp.call(['mkdir', run_name])
        sp.call(['mkdir', run_name + '/model_files'])
        print('Starting {}'.format(run_name))

        with open(run_name + '/samples.csv', 'w') as f:
            np.savetxt(f, (np.append([param[0] for param in to_vary], 'lnprob\n'),),
                delimiter=',', fmt='%s')
        pos = [[param[1] + param[2]*np.random.randn() for param in to_vary]
            for i in range(nwalkers)]

    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        print("Step {}".format(i))
        pos, chisum, blob = result
        with open(run_name + '/samples.csv', 'a') as f:
            np.savetxt(f, [np.append(pos[i], chisum[i]) for i in range(nwalkers)], delimiter=',')
        sp.call('rm -rf {}/model_data/*'.format(run_name), shell=True)

    pool.close()

#==============================================================================#
# Sandbox
#==============================================================================#

run_emcee(nsteps=10000, nwalkers=50, run_name='run6_50walkers_10params', to_vary = [
    ('m_disk',             -7.55,         0.05,       (-np.inf, np.inf)),
    ('sb_law',             2.3,           2,          (-5.,     10.)),
    ('scale_factor',       0.05,          0.05,       (0,       np.inf)),
    ('r_in',               8.8,           5,          (0,       np.inf)),
    ('d_r',                31.5,          10,         (0,       np.inf)),
    ('inc',                90,            1,          (0,       np.inf)),
    ('pa',                 128.48,        0.1,        (0,       360)),
    ('mar_starflux',       3.67e-4,       1e-4,       (0,       np.inf)),
    ('aug_starflux',       1.23e-4,       1e-4,       (0,       np.inf)),
    ('jun_starflux',       2.62e-4,       1e-4,       (0,       np.inf))])

# Display clean images and residuals for each observation
# standard_model = Model(params.values(), band6_observations)
# for obs in standard_model.observations:
#     standard_model.clean(obs)
#     raw_input('cool?')
#     standard_model.residuals(obs)
#     raw_input('cool?')
