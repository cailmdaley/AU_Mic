#This code is my attempt to write my own Metropolis-Hastings algorithm, and apply it to fake data.
import numpy as np
import math
import numpy.random as rand
import matplotlib.pyplot as plt
import scipy.optimize as op
from matplotlib.ticker import MaxNLocator
import emcee
import corner

#Create randomly sampled parabola
x = 10 * rand.random((100,)) - 5
gaussian = rand.normal(0, 1.6, 100)
data = 3*x**2 + gaussian + 1


def chi(a,b):
    model = a*x**2 + b
    chi = np.sum(((data-model)**2)/1.6**2)
    return chi

#emcee
def lnlike(theta, x, y, yerr):
    a, b = theta
    model = a*x**2 + b
    inv_sigma2 = 1.0/(yerr**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2))

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [2.5, 0.75], args=(x, data, gaussian)) #result["x"]

def lnprior(theta):
    a, b = theta
    if -10. < a < 10. and -10. < b < 10.0:
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lnlike(theta, x, y, yerr) + lp


ndim, nwalkers = 2, 100
#pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
pos = [[1.5, 5.] + 5*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, data, 1.6))
sampler.run_mcmc(pos, 500)


#burn in plot
plt.close()
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(3, color="#888888", lw=2)
axes[0].set_ylabel("$a$")

axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(1, color="#888888", lw=2)
axes[1].set_ylabel("$b$")
axes[1].set_xlabel("step number")

fig.tight_layout(h_pad=0.0)
plt.show()


#corner plot
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$a$", "$b$"],truths=[3, 1, np.log(gaussian)])
plt.show()









#M-H Algorithim
#a = [5.]
#b= [6.]
#steps = 10000
#a_accept = []
#alength = []
#b_accept = []
#blength = []
#chisquared = []
#sigma_a = 0.05
#sigma_b = 0.5
#for i in range(steps):
#    if rand.random(1) < 0.5:
#        alength.append(1)
#        a.append(rand.normal(a[-1], sigma_a))
#        b.append(b[-1])
#        deltachi = chi(a[-1], b[-1]) - chi(a[-2], b[-2])
#        if deltachi > 0:
#            if rand.random(1) <  np.exp(-deltachi/2):
#                a[-1] = a[-1]
#                a_accept.append(1)
#            else: a[-1] = a[-2]
#        else:
#            a[-1] = a[-1]
#            a_accept.append(1)
#    else:
#        blength.append(1)
#        b.append(rand.normal(b[-1], sigma_b))
#        a.append(a[-1])
#        deltachi = chi(a[-1], b[-1]) - chi(a[-2], b[-2])
#        if deltachi > 0:
#            if rand.random(1) <  np.exp(-deltachi/2):
#                b[-1] = b[-1]
#                b_accept.append(1)
#            else: b[-1] = b[-2]
#        else:
#            b[-1] = b[-1]
#            b_accept.append(1)
#    chisquared.append(chi(a[-1], b[-1]))
#a_acceptance = len(a_accept)/np.float(len(alength))
#b_acceptance = len(b_accept)/np.float(len(blength))
#acceptance = (len(a_accept) + len(b_accept))/np.float(steps)
#print a_acceptance
#print b_acceptance
#print acceptance

#Model & Data Plotting
def fitplot():
    plt.close()
    model = a[-1]*x**2 + b[-1]
    xindices = x.argsort()
    xsorted = x[xindices]
    modelsorted = model[xindices]
    plt.scatter(x, data,)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot (xsorted, modelsorted)
    plt.show(False)

#Chi^2 Plotting
def chiplot():
    plt.close()
    plt.loglog(range(steps),chisquared)
    plt.xlabel('steps (log)')
    plt.ylabel(r'$\chi^2$ (log)')
    plt.show(False)

#Histogram
def histplot():
    plt.close()
    binsize = 20
    a_hist = np.histogram(a, bins=binsize)
    height = math.ceil(max(a_hist[1]) - math.floor(min(a_hist[1])))
    height = height/len(a_hist[0])

    ax1= plt.subplot(2,1,1)
    ax1.barh(a_hist[1][:-1],a_hist[0], height=height)
    plt.xlabel('a-value')
    plt.ylabel('frequency')
    ax2=plt.subplot(2,1,2)
    ax2.hist(b,bins=binsize)
    plt.xlabel('b-value')
    plt.ylabel('frequency')
    plt.show(False)
