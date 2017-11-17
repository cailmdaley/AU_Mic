import pandas as pd
from astropy.io import fits
import scipy.stats
import numpy as np
from glob import glob

# create df
run_info = pd.DataFrame(columns=['samples', 'k',  'chi', 'red_chi', 'AICc', 'AIC_sigma', 'd_BIC'], dtype=float)
run_info.name = 'run'

# get # of visibilities
uvfs = glob('subtracted_obs_files/*.uvf')
nvis = np.sum([fits.open(uvf)[0].data.shape[0] for uvf in uvfs])


# get get run directories and names
runs = []
for i in range(6,100):
    try: runs.append(glob('run{}'.format(i))[0])
    except IndexError: pass
run_names = [glob(run + '_*')[0][:-3] for run in runs]

# pull nsamples, number of free params, and min chi^2 from each chain
for run, run_name in zip(runs, run_names):
    chain = pd.read_csv('{}/{}_chain.csv'.format(run,run))

    run_info.loc[run_name, 'samples'] = chain.shape[0]
    run_info.loc[run_name, 'k'] = chain.shape[1]
    run_info.loc[run_name, 'chi'] = chain.lnprob.max() * -2
run_info.red_chi = run_info.chi / nvis

# calculate AIC --> AICc
AIC = 2 * run_info.k + run_info.chi
run_info.AICc = AIC + (2*run_info.k * (run_info.k + 1)) / (nvis - run_info.k - 1)
# AICc differs from AIC by only 0.001, but what the hell let's use it.

# calculate model relative liklehood, converto to guassian sigma
rel_like = np.exp((run_info.AICc['run19_fiducial'] - run_info.AICc)/2.)
gauss = scipy.stats.norm()
for row in run_info.index:
    try:
        run_info.loc[row, 'AIC_sigma'] =  gauss.interval(1. - rel_like[row])[1]
    except ValueError:
        run_info.loc[row, 'AIC_sigma'] = np.nan

# calculate delta BIC
run_info.d_BIC = np.log(nvis)*run_info.k + run_info.chi
run_info.d_BIC -= run_info.d_BIC['run19_fiducial']

print(run_info.sort_values('AIC_sigma', ascending=True).round(3))
