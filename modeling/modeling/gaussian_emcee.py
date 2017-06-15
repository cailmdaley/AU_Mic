import numpy.random as rand
import numpy as np
import pandas as pd
import seaborn as sns
import emcee

    
xs = 10 * rand.random((100,)) - 5
xs.sort()
noise = np.array([1.6 for x in xs])
gaussian = 1 * rand.normal(0, noise)
ys = 15*np.exp(-30*xs**2) + 30*np.exp(-30*(xs-2)**2) + gaussian
data = (xs, ys, noise)

def model_function(x, (A1, A2, width, x0, x1)): 
    y = A1*np.exp(-width*(x-x0)**2) + A2 * np.exp(-width*(x-x1)**2)
    return y

def lnlike(theta, data):
    x, y, yerr = data
    model = model_function(x, theta)
    return -0.5 * (np.sum((y-model)**2 / yerr**2))
    
def lnprior(theta):
    A1, A2, width, x0, x1 = theta
    
    if A1 >=0 and A2 >= 0 and width >= 0 and abs(x0) < 10 and abs(x1) < 10:
        return 0.0
    return -np.inf

def lnprob(theta, data):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, data)
    
def fitplot(sampler):
    """Plot model fit on top of data"""
    plt.close()

    # Plot data
    x, y, yerr = sampler.args[0]
    plt.errorbar(x, y, yerr=yerr, fmt='.', elinewidth=0.5, capsize=0.5)

    # Plot random samples from posterior
    for sample in sampler.flatchain[rand.randint(len(sampler.flatchain), size=100)]:
        plt.plot(x, model_function(x, sample), 'k-', lw=0.5, alpha=0.7)

    # Plot final fit
    model = model_function(x, sampler.flatchain.mean(axis=0))
    plt.plot(x, model, 'r', lw=2, label='random posterior samples')
    plt.xlabel('$t$')
    plt.ylabel('$I$')
    plt.title("Data Fit")

    sns.despine()
    plt.show(False)
    # plt.savefig("MH_lightcurve1_fitplot.png")
    
def pairplot(sampler, param_names):
    """ Plot 'corner plot' of fit"""
    plt.close()
    posterior = pd.DataFrame(sampler.flatchain, columns=param_names)

    # cmap = sns.cubehelix_palette(as_cmap=True, start=2.3, dark=0, light=1, reverse=True)
    cmap = "Blues"
    corner = sns.PairGrid(posterior, diag_sharey=False, despine=False)
    corner.map_diag(sns.kdeplot)
    corner.map_lower(sns.kdeplot, cmap=cmap, n_levels=5, shade=True)
    corner.map_upper(plt.scatter, s=0.3)
    

    plt.subplots_adjust(top=0.9)
    corner.fig.suptitle("Corner Plot")
    plt.show(False)
    # plt.savefig('MH_lightcurve1_pairgrid.png')
    
    
ndim, nwalkers = 5, 10
pos = [[14.80363748,  30.51691964, 31.97330834, 0.04468592, 2.03675352] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(data,))
sampler.run_mcmc(pos, 500)
sampler.acceptance_fraction
sampler.flatchain.mean(axis=0)

fitplot(sampler)    
pairplot(sampler, ('A1', 'A2', 'width', 'x0', 'x1'))

# np.shape(*sampler.args)
# plt.scatter(xs,ys)    
# plt.scatter(xs, model_function(xs,sampler.chain.mean(axis=(0,1))))
# plt.show()
