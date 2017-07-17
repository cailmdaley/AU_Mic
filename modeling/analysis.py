import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

run_name = 'run4-16walkers_7params'
posterior = pd.read_csv(run_name + '.csv')
posterior
    
def walker_evolution(run_name, nwalkers):
    plt.close()
    
    posterior = pd.read_csv(run_name + '.csv')
    posterior.index /= nwalkers
    
    axes = posterior.iloc[0::nwalkers].plot(figsize=(7, 2*(len(posterior.columns)-1)), 
        subplots=True, title=run_name + ' evolution',
        color='black', alpha=0.5) 
    for i in range(nwalkers-1):
        posterior.iloc[i+1::nwalkers].plot(subplots=True, ax=axes, legend=False, 
            color='black', alpha=0.5)
        
    #quickly find misbehaving walkers
    # posterior.iloc[0::nwalkers].plot(subplots=True, ax=axes, color='red', alpha=0.5, legend=False)
    
    plt.savefig('{}.evolution.png'.format(run_name))
    plt.show()
    
def corner(run_name, nwalkers, text_specs, burn_in=0, bad_walkers=[]):
    """ Plot 'corner plot' of fit"""
    plt.close()
    
    # read in and groom posterior: remove burn-in, lnprob column, and bad walkers
    samples = pd.read_csv(run_name + '.csv')
    posterior = samples.iloc[burn_in*nwalkers:]
    posterior = posterior.drop([row for row in posterior.index if row%nwalkers in bad_walkers])
    posterior=posterior.drop(posterior.columns[-1], axis=1)
    posterior.rename(columns={'m_disk' : '$\log$ m_disk'}, inplace=True)

    # make corner plot
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=True)
    corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
    corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=3, shade=False)
    
    corner.fig.text(*text_specs[:-1], s=posterior.mean(), fontsize = text_specs[-1])
    
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
    corner.fig.suptitle(title, fontsize=25)
    corner.savefig(title + '.png')
    plt.show(False)
    
    # posterior.hist(bins=20)
    # plt.savefig('{}.hist.png'.format(run_name)); plt.show()
    
    
# walker_evolution('run3_8walkers_3params', nwalkers=8)
# corner('run3_8walkers_3params', nwalkers=8
#     text_specs = (0.45,0.7, 15), burn_in=50, bad_walkers=[7])

walker_evolution('run4-16walkers_7params', nwalkers=16)
corner('run4-16walkers_7params', nwalkers=16, text_specs = (0.45,0.7, 15),
    burn_in=600)
