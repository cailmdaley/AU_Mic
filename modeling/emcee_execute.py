import numpy as np
import emcee_fixed_starflux
from astrocail import fitting

fitting.run_emcee(run_name='run7', nsteps=10000, nwalkers=18, 
    lnprob = emcee_fixed_starflux.lnprob, to_vary = [
    ('m_disk',             -7.55,         0.05,       (-np.inf, np.inf)),
    ('sb_law',             2.3,           2,          (-5.,     10.)), 
    ('scale_factor',       0.05,          0.03,       (0,       np.inf)),
    ('r_in',               8.8,           5,          (0,       np.inf)),
    ('d_r',                31.5,          10,         (0,       np.inf)),
    ('inc',                90,            1,          (0,       np.inf)),
    ('pa',                 128.48,        0.1,        (0,       360)),
    ('starflux',           2.50e-4,       1e-4,      (0,       np.inf))])
