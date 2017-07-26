from emcee.utils import MPIPool
from aumic_obs_and_model import *

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
    
    # create model and get chi^2
    model = Model(params.values(), observations=band6_observations)
    model.chis=[]
    for date, starflux in zip(model.observations, model.starfluxes):
        for obs in date:
            model.obs_sample(obs, starflux)
            model.get_chi(obs)
    return -0.5 * sum(model.chis)
    
def run_mcmc(nsteps, nwalkers, run_name, to_vary):
    
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    # run sampler chain
    ndim = len(to_vary)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(run_name[:4], to_vary), pool=pool)
    
    try:
        df = pd.read_csv(run_name + '.csv')
        pos = np.array(df.iloc[-nwalkers:, :-1])
    except IOError:
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

import time
start=time.time()
run_mcmc(nsteps=10000, nwalkers=26, run_name='run5_26walkers_10params', to_vary = [
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
print('Run completed in {} hours'.format((time.time()-start) / 3600.))

# Display clean images and residuals for each observation
# standard_model = Model(params.values(), band6_observations)
# for obs in standard_model.observations:
#     standard_model.clean(obs)
#     raw_input('cool?')
#     standard_model.residuals(obs)
#     raw_input('cool?')
