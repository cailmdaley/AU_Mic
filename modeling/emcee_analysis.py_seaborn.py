# from model_plotting import FitsImage
# from aumic_obs_and_model import *
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

def rename(samples):
    samples[['mar_starflux', 'jun_starflux', 'aug_starflux']] *= 1e6
    samples['d_r'] += samples['r_in']
    samples.rename(columns={
        'm_disk'       : '$\log$ m_disk',
        'r_in'         : 'r$_{in}$',
        'd_r'          : 'r$_{out}$',
        'mar_starflux' : 'Mar starflux ($\mu$Jy)',
        'aug_starflux' : 'Aug starflux ($\mu$Jy)',
        'jun_starflux' : 'Jun starflux ($\mu$Jy)'}, inplace=True)

def walker_evolution(run_name, nwalkers):
    plt.close()

    samples = pd.read_csv(run_name + '.csv')
    rename(samples)
    last_step = samples.iloc[-nwalkers:]
    bad_walkers = last_step[last_step['lnprob'] == -np.inf].index % nwalkers
    samples['step'] = samples.index // nwalkers

    color = 'red' if 0 in bad_walkers else 'black'
    axes = samples.iloc[0::nwalkers].plot(x='step',
        figsize=(7, 2*(len(samples.columns)-1)), subplots=True,
        title=run_name + ' evolution', color=color, alpha=0.5)
    for i in range(nwalkers-1):
        color = 'red' if i+1 in bad_walkers else 'black'
        samples.iloc[i+1::nwalkers].plot(x='step',
            subplots=True, ax=axes, legend=False, color=color, alpha=0.5)

    # remove bad walkers before taking mean
    samples.drop([row for row in samples.index if row%nwalkers in bad_walkers], inplace=True)
    samples.set_index('step', inplace=True)

    walker_means = pd.DataFrame([samples.loc[i].mean() for i in range(samples.index[-1])])
    walker_means.plot(subplots=True, ax=axes, legend=False, color='forestgreen', ls='--')

    plt.savefig('{}.evolution.png'.format(run_name), dpi=700)
    # plt.show()


def corner_plot(run_name, nwalkers, stat_specs, burn_in=0, bad_walkers=[]):
    """ Plot 'corner plot' of fit"""
    plt.close()

    # read in samples
    samples = pd.read_csv(run_name + '.csv')
    rename(samples)

    # cut out burn in and bad walkers
    # posterior = samples.iloc[-100:, :-1]
    posterior = samples.iloc[burn_in*nwalkers:, :-1]

    last_step = samples.iloc[-nwalkers:]
    bad_walkers = last_step[last_step['lnprob'] == -np.inf].index % nwalkers
    posterior.drop([row for row in posterior.index if row%nwalkers in bad_walkers], inplace=True)
    print('walkers {} removed from posterior.'.format(tuple(bad_walkers)))

    # make corner plot
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=True)
    corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=False)

    # get best_fit and posterior statistics
    stats = posterior.describe().drop(['count', 'min', 'max'])
    stats.loc['best fit'] = samples.drop('lnprob', 1).loc[samples['lnprob'].idxmax()]
    stats = stats.iloc[[-1]].append(stats.iloc[:-1])
    print(stats.round(2).to_latex())

    table_ax = corner.fig.add_axes([0,0,1,1], frameon=False)
    table_ax.axis('off')
    left, bottom = stat_specs
    pd.plotting.table(table_ax, stats.round(2), bbox=[left, bottom, 1-left, .12], edges='open', colLoc='right')

    # hide upper triangle, so that it's a conventional corner plot
    for i, j in zip(*np.triu_indices_from(corner.axes, 1)):
        corner.axes[i, j].set_visible(False)

    for ax in corner.axes.flat:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))

    if burn_in == 0:
        title = run_name + '.corner_ungroomed'
    else:
        title = run_name + '.corner_groomed'

    plt.subplots_adjust(top=0.9)
    corner.fig.suptitle(r'{} Parameters, {} Walkers, {} Steps $\to$ {} Samples'
        .format(posterior.shape[1], nwalkers, posterior.shape[0]//nwalkers, posterior.shape[0], fontsize=25))
    corner.savefig(title + '.png', dpi=700)

# def view_best_model(run_name):
#
#     samples = pd.read_csv(run_name + '.csv')
#     best_fit_params = samples.drop('lnprob', 1).loc[samples['lnprob'].idxmax()]
#
#     for param in best_fit_params.index:
#         params[param] = best_fit_params[param]
#     params['m_disk'] = 10**params['m_disk']
#     params['d_r'] += params['r_in']
#
#     # fig = plt.subplots(1, 2, sharex=False, sharey=False,
#     #     figsize=(11.6, 6.5))[0]
#     # plt.subplots_adjust(wspace=-0.0)
#
#     model = Model(params=params.values(), observations=band6_observations,
#         name='best_fit')
#
#     for spw in band6_observations[2]:
#         model.obs_sample(spw, model.starfluxes[2])
#
#     cat_files = ','.join(['model_data/'+model.name+obs.name+'.vis' for obs in band6_observations[2]])
#     sp.call(['uvcat', 'vis={}'.format(cat_files), 'out=model_data/best_fit.vis'])
#
#     # Clean down to half the observation rms
#     sp.call(['invert',
#         'vis=model_data/{}.vis'.format('best_fit'),
#         'map=model_data/{}.mp'.format('best_fit'),
#         'beam=model_data/{}.bm'.format('best_fit'),
#         'cell=0.03arcsec', 'imsize=512', 'options=systemp,mfs', 'robust=2'])
#     sp.call(['clean',
#         'map=model_data/{}.mp'.format('best_fit'),
#         'beam=model_data/{}.bm'.format('best_fit'),
#         'out=model_data/{}.cl'.format('best_fit'),
#         'niters=100000', 'cutoff={}'.format(obs.rms/2)])
#     sp.call(['restor',
#         'map=model_data/{}.mp'.format('best_fit'),
#         'beam=model_data/{}.bm'.format('best_fit'),
#         'model=model_data/{}.cl'.format('best_fit'),
#         'out=model_data/{}.cm'.format('best_fit')])
#     sp.call(['fits', 'op=xyout',
#         'in=model_data/{}.fits'.format('best_fit'),
#         'out=model_data/{}.im'.format('best_fit')], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
#
#
#
#     test = FitsImage('model_data/' + model.name + '.fits', rms=1.5e-05, text=[[4.6,4.0, 'Best Fit Model']],
#         cbspace=[1, 0.2])
#     test.im.max()

# all_natural = Observation('../../cleans/current/aumic_band6_all_natural.fits',
#     1.4494822607957758e-05, fig=fig, pos=(0, num),
#     text=[[4.6, 4.0, 'ALMA 1.4 mm']])
#         #   [4.6, 3.0, 'natural weighting']])
# view_best_model('run5_26walkers_10params')
# plt.savefig('best_fit_model.png', dpi=700)
# plt.show()

# walker_evolution('run3_8walkers_3params', nwalkers=8)
# corner('run3_8walkers_3params', nwalkers=8
#     text_specs = (0.45,0.7, 15), burn_in=50, bad_walkers=[7])

# walker_evolution('run4-16walkers_7params', nwalkers=16)
# corner_plot('run4-16walkers_7params', nwalkers=16, burn_in=600, text_specs=(0.19,0.8, 16))


# walker_evolution('run5_26walkers_10params', nwalkers=26)
# corner_plot('run5_26walkers_10params', nwalkers=26, burn_in=400, stat_specs=(.18, .82))

#walker_evolution('run6_50walkers_10params', nwalkers=50)
#corner_plot('run6_50walkers_10params', nwalkers=50, burn_in=1000, stat_specs=(.18, .82))

walker_evolution('run7_constant_flux_18_walkers', nwalkers=18)
