import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

def pairplot(run_name):
    """ Plot 'corner plot' of fit"""
    plt.close()
    
    # run_name = 'run2-priors'
    # run_name = 'run1'
    posterior = pd.read_csv(run_name + '.csv')
    posterior['iteration'] = pd.cut(posterior.index.values, 2)
    # posterior.ix[1000:,'burn_in'] = 0

    
    # try:
    #     posterior['m_disk'] = 3.67 * 10**posterior['m_disk']
    # except KeyError:
    #     pass

    # cmap = sns.cubehelix_palette(as_cmap=True, start=2.3, dark=0, light=1, reverse=True)
    cmap = "Blues"
    corner = sns.PairGrid(posterior, hue='iteration', diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    corner.map_lower(sns.kdeplot, cut=0, cmap=cmap, n_levels=10, shade=True)
    corner.map_upper(plt.scatter)
    corner.add_legend()
    for ax in corner.axes.flat:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
    plt.subplots_adjust(top=0.95)
    corner.fig.suptitle("{}".format(run_name))
    
    corner.savefig('{}.pairgrid.png'.format(run_name)); 
    plt.show(False)
    
    # posterior.hist(bins=20)
    # plt.savefig('{}.hist.png'.format(run_name)); plt.show()
    
    return posterior

test=pairplot('run2-priors')
# test=pairplot('run1')
