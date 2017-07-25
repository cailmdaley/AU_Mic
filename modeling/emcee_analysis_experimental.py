import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from sklearn.neighbors.kde import KernelDensity
from sklearn.model_selection import GridSearchCV
import pandas as pd
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.set_context("talk")

run_name = 'run5_26walkers_10params'
posterior = pd.read_csv(run_name + '.csv')


multi = posterior[['m_disk', 'sb_law']]

sns.kdeplot(multi)
plt.show()

bw_search = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(0, multi.std().mean()/10, 5)}, cv=20)
bw_search.fit(multi)
multi.shape[0]**(-1./(multi.shape[1]+4))
multi.std()
**()

xx, yy = np.meshgrid(np.linspace(*multi.iloc[:,0].quantile([0,1]), num=100), 
                     np.linspace(*multi.iloc[:,1].quantile([0,1]), num=100))
test = np.array([xx.ravel(), yy.ravel()]).T
kde=KernelDensity(bandwidth=multi.std().min()/10)
kde.fit(multi)
pdf = np.exp(kde.score_samples(multi)).reshape(len(xx), -1)

plt.contour(xx, yy, pdf )
plt.savefig('test.png')
plt.show()



def my_kde(df):
    mdisk = posterior['sb_law'].values.reshape(-1,1)
    bw_search = GridSearchCV(KernelDensity(), {'bandwidth': np.linspace(0, mdisk.std(), 5)})
    bw_search.fit(mdisk)
    kde = bw_search.best_estimator_.fit((mdisk))
    x = np.linspace(mdisk.min(), mdisk.max(), 1000)[:,None]
    pdf = np.exp(kde.score_samples(x))
    
    plt.plot(x, pdf)
    plt.show()
    
def walker_evolution(run_name, nwalkers):
    plt.close()
    
    samples = pd.read_csv(run_name + '.csv')
    last_step = samples.iloc[-nwalkers:]
    bad_walkers = last_step[last_step['lnprob'] == -np.inf].index % nwalkers
    samples.index //= nwalkers
    
    color = 'red' if 0 in bad_walkers else 'black'
    axes = samples.iloc[0::nwalkers].plot(
        figsize=(7, 2*(len(samples.columns)-1)), subplots=True, 
        title=run_name + ' evolution', color=color, alpha=0.5) 
    for i in range(nwalkers-1):
        color = 'red' if i+1 in bad_walkers else 'black'
        samples.iloc[i+1::nwalkers].plot(
            subplots=True, ax=axes, legend=False, color=color, alpha=0.5)
    walker_means = pd.DataFrame([samples.loc[i].mean() for i in range(samples.index[-1])])
    walker_means.plot(subplots=True, ax=axes, legend=False, color='forestgreen', ls='--')
        
    plt.savefig('{}.evolution.png'.format(run_name), dpi=700)
    plt.show()
    
    
def corner_plot(run_name, nwalkers, stat_specs, burn_in=0, bad_walkers=[]):
    """ Plot 'corner plot' of fit"""
    plt.close()
    
    # read in samples
    samples = pd.read_csv(run_name + '.csv')
    samples[['mar_starflux', 'jun_starflux', 'aug_starflux']] *= 1e6
    samples.rename(columns={
        'm_disk'       : '$\log$ m$_{\text{disk}}$',
        'd_r'          : '$\Delta$r',
        'mar_starflux' : 'Mar starflux ($\mu$Jy)',
        'aug_starflux' : 'Aug starflux ($\mu$Jy)',
        'jun_starflux' : 'Jun starflux ($\mu$Jy)'}, inplace=True)
    
    # cut out burn in and bad walkers
    # posterior = samples.iloc[burn_in*nwalkers:, :-1]
    posterior = samples.iloc[-1000:, :-1]
    
    last_step = samples.iloc[-nwalkers:]
    bad_walkers = last_step[last_step['lnprob'] == -np.inf].index % nwalkers
    posterior.drop([row for row in posterior.index if row%nwalkers in bad_walkers], inplace=True)
    print('walkers {} removed from posterior.'.format(tuple(bad_walkers)))
        
        
    # make corner plot
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot, cut=0)
    # hide upper triangle, so that it's a conventional corner plot
    for i, j in zip(*np.triu_indices_from(corner.axes, 1)):
        corner.axes[i, j].set_visible(False)
    for i, j in zip(*np.tril_indices_from(corner.axes, -1)):
        ax = corner.axes[i, j]
        plt.sca(ax)
        
        x_var = corner.x_vars[j]
        y_var = corner.y_vars[i]
        x_var
        y_var
        
        evaluate(scipy.stats.gaussian_kde, np.array([corner.data[x_var], corner.data[y_var]]))
        H, xedges,yedges=np.histogram2d(corner.data[x_var],corner.data[y_var],bins=20,normed=True)
        norm=H.sum() # Find the norm of the sum
        norm
        # Set contour levels
        contour1=0.99
        contour2=0.95
        contour3=0.68

        # Set target levels as percentage of norm
        target1 = norm*contour1
        target2 = norm*contour2
        target3 = norm*contour3

        # Take histogram bin membership as proportional to Likelihood
        # This is true when data comes from a Markovian process
        def objective(limit, target):
            w = np.where(H>limit)
            count = H[w]
            return count.sum() - target

        import scipy
        # Find levels by summing histogram to objective
        level1= scipy.optimize.bisect(objective, H.min(), H.max(), args=(target1,))
        level2= scipy.optimize.bisect(objective, H.min(), H.max(), args=(target2,))
        level3= scipy.optimize.bisect(objective, H.min(), H.max(), args=(target3,))
        level1

        # For nice contour shading with seaborn, define top level
        level4=H.max()
        levels=[level1,level2,level3,level4]
        
        sns.kdeplot(corner.data[x_var], corner.data[y_var], n_levels=levels)
        plt.show()
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3g'))
    plt.show()
    # corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=[.68, .95, .99], shade=True)
    # corner.map_lower(plt.scatter, s=1, color='#708090', alpha=0.2)
    # corner.map_lower(sns.kdeplot, cut=0, cmap='Blues', n_levels=[[.5,.8] for ax in corner.axes.flat], cbar=True, shade=False)
    
    # get best_fit and posterior statistics
    stats = posterior.describe().drop(['count', 'min', 'max'])
    stats.loc['best_fit'] = samples.drop('lnprob', 1).loc[samples['lnprob'].idxmax()]
    stats = stats.iloc[[-1]].append(stats.iloc[:-1])
    
    table_ax = corner.fig.add_axes([0,0,1,1], frameon=False)
    table_ax.axis('off')
    left, bottom = stat_specs
    pd.plotting.table(table_ax, stats.round(2), bbox=[left, bottom, 1-left, .12], colLoc='right')

        
                
    # for ax in corner.axes.flat:
        # plt.sca(ax  )
        # sns.kdeplot()    
    
    if burn_in == 0:
        title = run_name + '.corner_ungroomed'
    else:
        title = run_name + '.corner_groomed'
        
    plt.subplots_adjust(top=0.9)
    corner.fig.suptitle(r'{} Parameters, {} Walkers, {} Steps $\to$ {} Samples'
        .format(posterior.shape[1], nwalkers, posterior.shape[0]//nwalkers, posterior.shape[0], fontsize=25))
    corner.savefig(title + '.png')
    plt.show()
    
# walker_evolution('run3_8walkers_3params', nwalkers=8)
# corner('run3_8walkers_3params', nwalkers=8
#     text_specs = (0.45,0.7, 15), burn_in=50, bad_walkers=[7])

# walker_evolution('run4-16walkers_7params', nwalkers=16)
# corner_plot('run4-16walkers_7params', nwalkers=16, burn_in=600, text_specs=(0.19,0.8, 16))


# walker_evolution('run5_26walkers_10params', nwalkers=26)
corner_plot('run5_26walkers_10params', nwalkers=26, burn_in=400, stat_specs=(.18, .82))
