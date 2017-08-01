from aumic_obs_and_model import *
from emcee.utils import MPIPool
from aumic_obs_and_model import *

def make_fits(params, model):
    disk_params = params[:-1]
    PA = params[-1]

    model_disk = debris_disk.Disk(disk_params, obs=[300, 131, 300, 20])
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
def lnprob(theta, modelname, to_vary):
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
    model = Model(params.values(), observations=band6_observations, name=modelname)
    model.chis=[]
    for date, starflux in zip(model.observations, model.starfluxes):
        for obs in date:
            model.obs_sample(obs, starflux)
            model.get_chi(obs)
    return -0.5 * sum(model.chis)
def lnprior(theta):
    m_disk, sb_law, scale_factor, r_in, r_out, inc, pa = theta
    
def run_mcmc(nsteps, nwalkers, run_name, to_vary):
    
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    # run sampler chain
    ndim = len(to_vary)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(run_name[:4], to_vary), pool=pool)
    
    # try:
    #     posterior = pd.read_csv(run_name + '.csv')
    #     print('Resuming {} at step {}'.format(run_name, posterior.index[-1]//nwalkers))
    #     pos = np.array(posterior.iloc[-nwalkers:, :-1])
    #     with open(run_name + '.csv', 'a') as f:
    #         f.write('\n')
    # except IOError, IndexError:
    #     print('Starting {}'.format(run_name))
    #     with open(run_name + '.csv', 'w') as f:
    #         np.savetxt(f, (np.append([param[0] for param in to_vary], 'lnprob\n'),), 
    #             delimiter=',', fmt='%s')
    #     pos = [[param[1] + param[2]*np.random.randn() for param in to_vary] 
    #         for i in range(nwalkers)]
    try:
        samples = pd.read_csv(run_name + '.csv')
        print('Resuming {} at step {}'.format(run_name, samples.index[-1]//nwalkers))
        pos = np.array(samples.iloc[-nwalkers:, :-1])
        with open(run_name + '.csv', 'a') as f:
            f.write('\n')
    except IOError:
        print('Starting {}'.format(run_name))
        with open(run_name + '.csv', 'w') as f:
            np.savetxt(f, (np.append([param[0] for param in to_vary], 'lnprob\n'),), 
                delimiter=',', fmt='%s')
        pos = [[param[1] + param[2]*np.random.randn() for param in to_vary] 
            for i in range(nwalkers)] 
            
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        print("Step {}".format(i))
        pos, chisum, blob = result
        with open(run_name + '.csv', 'a') as f: 
            np.savetxt(f, [np.append(pos[i], chisum[i]) for i in range(nwalkers)], delimiter=',')
        sp.call('rm -rf model_data/*', shell=True)

    pool.close()
    
#==============================================================================#
# Sandbox
#==============================================================================#

run_mcmc(nsteps=10000, nwalkers=50, run_name='run6_50walkers_10params', to_vary = [
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
