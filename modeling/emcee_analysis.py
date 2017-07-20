import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

run_name = 'run5_26walkers_10params'
posterior = pd.read_csv(run_name + '.csv')
    
def walker_evolution(run_name, nwalkers):
    plt.close()
    
    samples = pd.read_csv(run_name + '.csv')
    bad_walkers = samples.iloc[-nwalkers:][samples['lnprob'] == -np.inf].index % nwalkers
    samples.index //= nwalkers
    
    axes = samples.iloc[0::nwalkers].plot(figsize=(7, 2*(len(samples.columns)-1)), 
        subplots=True, title=run_name + ' evolution',
        color='black', alpha=0.5) 
    for i in range(nwalkers-1):
        samples.iloc[i+1::nwalkers].plot(subplots=True, ax=axes, legend=False, 
            color='black', alpha=0.5)
        
    #quickly find misbehaving walkers
    for walker in bad_walkers:
        samples.iloc[walker::nwalkers].plot(subplots=True, ax=axes, color='red', alpha=0.5, legend=False)
    
    plt.savefig('{}.evolution.png'.format(run_name), dpi=700)
    plt.show()
    
def corner_plot(run_name, nwalkers, text_specs, burn_in=0, bad_walkers=[]):
    """ Plot 'corner plot' of fit"""
    plt.close()
    
    # read in samples
    samples = pd.read_csv(run_name + '.csv')
    samples[['mar_starflux', 'jun_starflux', 'aug_starflux']] *= 1e6
    samples.rename(columns={
        'm_disk'       : '$\log$ m_disk',
        'd_r'          : '$\Delta$r',
        'mar_starflux' : 'mar_starflux ($\mu$Jy)',
        'aug_starflux' : 'aug_starflux ($\mu$Jy)',
        'jun_starflux' : 'jun_starflux ($\mu$Jy)'}, inplace=True)
    
    # cut out burn in and bad walkers
    posterior = samples.iloc[burn_in*nwalkers:, :-1]
    blah=np.where(posterior.duplicated(keep=False)==True)[0]
    blah = np.where(posterior.duplicated()==True)[0] % nwalkers
    blah[:20]
    
    bad_walkers = samples.iloc[-nwalkers:][samples['lnprob'] == -np.inf].index % nwalkers
    posterior = posterior.drop([row for row in posterior.index if row%nwalkers in bad_walkers])
        
    # make corner plot
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=True)
    corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=False)
    
    # get best_fit and posterior statistics
    posterior_describe = posterior.describe().drop(['count', 'min', 'max'])
    best_fit = samples.iloc[samples['lnprob'].idxmax()].to_frame().transpose()
    best_fit.index = ['best_fit']
    stats = best_fit.append(posterior_describe)[best_fit.columns].fillna('')
    
    table_ax = corner.fig.add_axes([0,0,1,1], frameon=False)
    table_ax.axis('off')
    pd.plotting.table(table_ax, stats.round(2), bbox=[.24,.79,.76,.15], edges='open', 
        colLoc='right')

    # corner.fig.text(*text_specs[:-1],  fontsize = text_specs[-1]  ,
    #     s='Posterior statistics: {} steps, {} walkers, {} samples \n{}'
    #     .format(len(posterior)//nwalkers, nwalkers, len(posterior), stats_string), ha='left', multialignment='left')
    
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
    corner.savefig(title + '.png')
    
    
# walker_evolution('run3_8walkers_3params', nwalkers=8)
# corner('run3_8walkers_3params', nwalkers=8
#     text_specs = (0.45,0.7, 15), burn_in=50, bad_walkers=[7])

walker_evolution('run4-16walkers_7params', nwalkers=16)
corner_plot('run4-16walkers_7params', nwalkers=16, burn_in=600, text_specs=(0.19,0.8, 16))


walker_evolution('run5_26walkers_10params', nwalkers=26)
corner_plot('run5_26walkers_10params', nwalkers=26, text_specs=(0.19,0.8, 16))
