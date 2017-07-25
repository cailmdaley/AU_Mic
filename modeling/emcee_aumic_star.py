from aumic_obs_and_model import *
from emcee.utils import MPIPool

# define likelehood functions
# def lnprob(theta, modelname, to_vary):
#     """
#     For each parameter to be varied, return -infinity if the proposed value lies
#     outside the bounds for the parameter; if all parameters are ok, return 0.
#     """
#     
#     # reassign parameters to vary to those input by emcee
#     for i, free_param in enumerate(to_vary):
#         lower_bound, upper_bound = free_param[-1]
#         print(free_param)
#         print(theta[i])
#         
#         if lower_bound < theta[i] < upper_bound:
#             params[free_param[0]] = theta[i]
#         else: return -np.inf
#         
#     #logarithmic sampling
#     params['m_disk'] = 10**params['m_disk'] 
#     # since (k)evan's code takes r_out rather than d_r, add r_in to d_r to get r_out
#     params['d_r'] += params['r_in']
#     
#     # create model
#     model = Model(params.values(), observations=band6_observations, name=modelname)
#     
#     # return chi^2
#     print(-0.5 * sum(model.chis))
#     return -0.5 * sum(model.chis)

# define likelehood functions
def lnlike(theta, to_vary, observations):
    # default parameter values
    params = OrderedDict([
        ('temp_index',        -0.5),
        ('m_disk',            3.67e-08),
        ('sb_law',            2.3),
        ('r_in',              8.8),
        ('r_out',             40.3),
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
        ('mar_starflux',      0.000367),
        ('mar_starflux',      0.000367),
        ])
    
    # reassign parameters to vary to those input by emcee
    for i in range(len(to_vary)):
        params[to_vary[i][0]] = theta[i]
        
    params['m_disk'] = 10**params['m_disk']
    
    # create model
    model = Model(params.values(), observations=observations)
    
    # return chi^2
    return -0.5 * sum(model.chis)
def lnprior(theta):
    m_disk, sb_law, scale_factor, r_in, r_out, inc, pa = theta
    
    if -11. < m_disk and \
    -4. < sb_law < 10. and \
    0. < scale_factor < 2. and  \
    0. < r_in < r_out and \
    inc <= 90: 
        return 0.0
    else:
        return -np.inf
def lnprob(theta, to_vary, observations):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, to_vary, observations)
    
    
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
