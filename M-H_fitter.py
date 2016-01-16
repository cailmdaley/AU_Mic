#This code is my attempt to write my own Metropolis-Hastings algorithm, and apply it to fake data.
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
plt.close()

#Create randomly sampled parabola
x = 10 * rand.random((100,)) - 5
gaussian = 8*rand.normal(0, .2, 100)
data = 3*x**2 + gaussian + 1

def chi(a,b):
    model = a*x**2 + b
    chi = np.sum(((data-model)**2)/model)
    return chi

#M-H Algorithim
a = [5.]
b= [6.]
steps=1000
accept= []
chisquared = []
for i in range(steps):
    if rand.random(1) < 0.5:
        a.append(rand.normal(a[-1], 1))
        deltachi = chi(a[-2], b[-1]) - chi(a[-1], b[-1])
        if deltachi < 0:
            if rand.random(1) <  np.exp(-deltachi**2/2):
                a[-1] = a[-1]
                if i<30:
                    accept.append(1)
            else: a[-1] = a[-2]
        else:
            a[-1] = a[-1]
            if i<30:
                accept.append(1)
    else:
        b.append(rand.normal(b[-1], 1))
        deltachi = chi(a[-1], b[-2]) - chi(a[-1], b[-1])
        if deltachi < 0:
            if rand.random(1) <  np.exp(-deltachi**2/2):
                b[-1] = b[-1]
                if i<30:
                    accept.append(1)
            else: b[-1] = b[-2]
        else:
            b[-1] = b[-1]
            if i<30:
                accept.append(1)
    chisquared.append(chi(a[-1], b[-1]))
acceptance = len(accept)/30.
print acceptance

#Model & Data Plotting
model = a[-1]*x**2 + b[-1]
xindices = x.argsort()
xsorted = x[xindices]
modelsorted = model[xindices]
plt.scatter(x, data,)
plt.plot (xsorted, modelsorted)
plt.show(False)

#Chi^2 Plotting
# plt.plot(range(steps),chisquared)
# plt.show(False)
