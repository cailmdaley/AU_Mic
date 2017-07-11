import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

def pairplot(run_name):
    """ Plot 'corner plot' of fit"""
    posterior = pd.read_csv(run_name + '.csv')

    # cmap = sns.cubehelix_palette(as_cmap=True, start=2.3, dark=0, light=1, reverse=True)
    cmap = "Blues"
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot)
    corner.map_lower(sns.kdeplot, cmap=cmap, n_levels=5, shade=True)
    corner.map_upper(plt.scatter, s=0.7)
    
    plt.show(False)
    plt.savefig('pairgrid.png')
    
    return posterior

test=pairplot('test')
test.hist(bins=40)
plt.show()
