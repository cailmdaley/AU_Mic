import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

def walker_evolution(run_name, nwalkers):
    plt.close()
    run_name = 'run3'
    posterior = pd.read_csv(run_name + '.csv')
    posterior.index /= nwalkers
    
    # posterior.reset_index(inplace=True)
    # posterior.rename(columns={'index' : 'walker'}, inplace=True)
    # posterior
    # posterior['walker'] %= nwalkers
    
    axes = posterior.iloc[0::nwalkers].plot(subplots=True, color='black', alpha=0.5)
    for i in range(nwalkers):
        posterior.iloc[i+1::nwalkers].plot(subplots=True, ax=axes, legend=False, color='black', alpha=0.5)
        
    #quickly find misbehaving walkers
    # posterior.iloc[7::nwalkers].plot(subplots=True, ax=axes, color='red', alpha=0.5, legend=False)
    
    plt.show()
         
def corner(run_name, nwalkers, burn_in=0, bad_walkers=[]):
    """ Plot 'corner plot' of fit"""
    plt.close()
    
    # run_name = 'run3'
    samples = pd.read_csv(run_name + '.csv' )
    posterior = samples.iloc[burn_in*nwalkers:]
        
    posterior = posterior.drop([row for row in posterior.index if row%nwalkers in bad_walkers])
    posterior=posterior.drop('lnprob', axis=1)
    posterior.rename(columns={'m_disk' : '$\log$ m_disk'}, inplace=True)
    print(posterior.mean())

    cmap = "Blues"
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    corner.map_lower(sns.kdeplot, cut=0, cmap=cmap, n_levels=3, shade=True)
    corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
    corner.map_lower(sns.kdeplot, cut=0, cmap=cmap, n_levels=3, shade=False)
    
    for i, j in zip(*np.triu_indices_from(corner.axes, 1)):
        corner.axes[i, j].set_visible(False)
    
    for ax in corner.axes.flat:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    plt.subplots_adjust(top=0.95)
    corner.fig.suptitle("{}".format(run_name))
    
    corner.savefig('{}.pairgrid.png'.format(run_name), dpi=700); 
    plt.show(False)
    
    # posterior.hist(bins=20)
    # plt.savefig('{}.hist.png'.format(run_name)); plt.show()
    
    return posterior
    
posterior.mean()
test=corner('run3', nwalkers=8, burn_in=50, bad_walkers=[7])
# walker_evolution('run3', nwalkers=8)
# test=pairplot('run1')
