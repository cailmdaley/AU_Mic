# This code is my attempt to write my own Metropolis-Hastings algorithm,
# and apply it to fake data.
import numpy as np
import math
import numpy.random as rand
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
import scipy.optimize as op
from matplotlib.ticker import MaxNLocator
# import emcee
# import corner


def chi(data, f, params):
    model = f(data[0], params)
    chi_squared = np.sum(((data[1] - model)**2) / data[2]**2)
    return chi_squared

class fit_info:
    def __init__(self, data, model, params, guesses, steps, chis):
        self.data = data
        self.model = model
        self.params = params
        self.guesses = guesses
        self.steps = steps
        self.chis = chis

# M-H Algorithim
def mh_fit(data, model_function, prior, sigmas, steps=10000):
    tries = np.zeros(len(prior))
    accepts = np.zeros(len(prior))
    chis = []
    guesses = [prior]

    sigma_a = 0.05
    sigma_b = 0.5
    
    for i in range(steps):
        param = rand.randint(len(prior))
        tries[param] += 1

        prev_guess = guesses[-1][:]
        current_guess = guesses[-1][:]
        current_guess[param] = rand.normal(prev_guess[param], sigmas[param])

        deltachi = chi(data, model_function,
                       current_guess) - chi(data, model_function, prev_guess)
                       
        if deltachi > 0 and rand.random(1) > np.exp(-deltachi / 2):
            guesses.append(prev_guess)
        else:
            guesses.append(current_guess)    
            accepts[param] += 1

        chis.append(chi(data, model_function, guesses[-1]))
        
    params = guesses[-1]
    model = model_function(data[0], params)    
    acceptances = accepts/tries
    tot_acceptance = np.sum(accepts)/np.sum(tries)
    
    print("Paremeter acceptance rate: {}".format(acceptances))
    print("Total acceptance rate: {}".format(tot_acceptance))
    print("Final chi^2: {}".format(chis[-1]))


    return fit_info(data=data, model=model, params=params, guesses=guesses, steps=steps, chis=chis)

# Model & Data Plotting
def fitplot(fit):
    plt.close()
    xindices = fit.data[0].argsort()
    xsorted = fit.data[0][xindices]
    modelsorted = fit.model[xindices]
    plt.scatter(fit.data[0], fit.data[1])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(xsorted, modelsorted)
    plt.show(False)

# Chi^2 Plotting


def chiplot(fit):
    plt.close()
    plt.loglog(range(fit.steps), fit.chis)
    plt.xlabel('steps (log)')
    plt.ylabel(r'$\chi^2$ (log)')
    plt.show(False)

# Histogram


def histplot(fit):
    plt.close()
    
    sns.jointplot(np.array([guess[0] for guess in fit.guesses]), np.array([guess[1] for guess in fit.guesses]), kind='resid')
    plt.xlabel('a-value')
    plt.ylabel('b-value')
    plt.show(False)
    sns.jointplot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create randomly sampled parabola
xs = 10 * rand.random((100,)) - 5
noise = np.array([1.6 for x in xs])
gaussian = 1 * rand.normal(0, noise, 100)
data = 3 * xs**2 + gaussian + 1
data = (xs, data, noise)
def model_function(x, (a, b)): return a * x**2 + b


fit = mh_fit(data, model_function, [2, 2], (0.05, 0.5))

fitplot(fit)
chiplot(fit)
histplot(fit)
