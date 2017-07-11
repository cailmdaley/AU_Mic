import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

run_name = 'run1'
def pairplot(run_name):
    """ Plot 'corner plot' of fit"""
    plt.close()
    
    posterior = pd.read_csv(run_name + '.csv')
    
    # try:
    #     posterior['m_disk'] = 3.67 * 10**posterior['m_disk']
    # except KeyError:
    #     pass

    # cmap = sns.cubehelix_palette(as_cmap=True, start=2.3, dark=0, light=1, reverse=True)
    cmap = "Blues"
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot)
    corner.map_lower(sns.kdeplot, cmap=cmap, n_levels=10, shade=True)
    corner.map_upper(plt.scatter, s=1.5)
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
test=pairplot('run1')
